#include "src/models/Acid/Acid2dSolver.hpp"

#include "adolc/sparse/sparsedrivers.h"
#include "adolc/drivers/drivers.h"
#include <iomanip>

using namespace acid2d;
using std::vector;
using std::ofstream;
using std::endl;
using std::map;
using std::setprecision;

Acid2dSolver::Acid2dSolver(Acid2d* _model) : AbstractSolver<Model>(_model)
{
	y = new double[var_size * size];

	const size_t strNum = var_size * model->cellsNum;
	ind_i = new int[mesh::stencil * var_size * strNum];
	ind_j = new int[mesh::stencil * var_size * strNum];
	cols = new int[strNum];
	a = new double[mesh::stencil * var_size * strNum];
	ind_rhs = new int[strNum];
	rhs = new double[strNum];

	options[0] = 0;          /* sparsity pattern by index domains (default) */
	options[1] = 0;          /*                         safe mode (default) */
	options[2] = 0;          /*              not required if options[0] = 0 */
	options[3] = 0;          /*                column compression (default) */

	P.open("snaps/P.dat", ofstream::out);
	S.open("snaps/S.dat", ofstream::out);
	qcells.open("snaps/Q.dat", ofstream::out);

	CHOP_MULT = 0.1;
	MAX_SAT_CHANGE = 1.0;
	
	CONV_W2 = 1.e-4;		CONV_VAR = 1.e-10;
	MAX_ITER = 20;
}
Acid2dSolver::~Acid2dSolver()
{
	delete[] y;
	delete[] ind_i, ind_j, ind_rhs;
	delete[] cols;
	delete[] a, rhs;

	P.close();
	S.close();
	qcells.close();
}
void Acid2dSolver::writeData()
{
	double p = 0.0, s = 0.0, q = 0.0;

	qcells << cur_t * t_dim / 3600.0;

	map<size_t, double>::iterator it;
	for (it = model->Qcell.begin(); it != model->Qcell.end(); ++it)
	{
		const auto data = (*model)[it->first].u_next;
		p += data.p * model->P_dim;
		s += data.s;
		if (model->leftBoundIsRate)
			qcells << "\t" << it->second * model->Q_dim * 86400.0;
		else
		{
			qcells << "\t" << model->getRate(it->first) * model->Q_dim * 86400.0;
			q += model->getRate(it->first);
		}
	}
	P << cur_t * t_dim / 3600.0 << "\t" << p / (double)(model->Qcell.size()) << endl;
	S << cur_t * t_dim / 3600.0 << "\t" << s / (double)(model->Qcell.size()) << endl;

	if (model->leftBoundIsRate)
		qcells << "\t" << model->Q_sum * model->Q_dim * 86400.0 << endl;
	else
		qcells << "\t" << q * model->Q_dim * 86400.0 << endl;
}
void Acid2dSolver::control()
{
	writeData();

	if (cur_t >= model->period[curTimePeriod])
	{
		curTimePeriod++;
		model->ht = model->ht_min;
		model->setPeriod(curTimePeriod);
	}

	if (model->ht <= model->ht_max && iterations < 6)
		model->ht = model->ht * 1.5;
	else if (iterations > 6 && model->ht > model->ht_min)
		model->ht = model->ht / 1.5;

	if (cur_t + model->ht > model->period[curTimePeriod])
		model->ht = model->period[curTimePeriod] - cur_t;

	cur_t += model->ht;
}
void Acid2dSolver::start()
{
	int counter = 0;
	iterations = 8;

	fillIndices();
	solver.Init(var_size * model->cellsNum, 1.e-15, 1.e-15);

	model->setPeriod(curTimePeriod);

	while (cur_t < Tt)
	{
		control();
		model->snapshot_all(counter++);
		doNextStep();
		copyTimeLayer();
		cout << "---------------------NEW TIME STEP---------------------" << endl;
		cout << setprecision(6);
		cout << "time = " << cur_t << endl;
	}
	model->snapshot_all(counter++);
	writeData();
}
void Acid2dSolver::copySolution(const paralution::LocalVector<double>& sol)
{
	for (size_t i = 0; i < size; i++)
	{
		auto& var = (*model)[i].u_next;
		var.m += sol[i * var_size];
		var.p += sol[i * var_size + 1];
		var.s += sol[i * var_size + 2];
		var.xa += sol[i * var_size + 3];
		var.xw += sol[i * var_size + 4];
	}
}
void Acid2dSolver::checkStability()
{
	auto barelyMobilLeft = [this](double s_cur, double s_crit) -> double
	{
		return s_crit + fabs(s_cur - s_crit) * CHOP_MULT;
	};
	auto barelyMobilRight = [this](double s_cur, double s_crit) -> double
	{
		return s_crit - fabs(s_crit - s_cur) * CHOP_MULT;
	};
	auto checkCritPoints = [=, this](auto& next, auto& iter, auto& props)
	{
		// Oil
		if ((1.0 - next.s - props.s_oc) * (1.0 - iter.s - props.s_oc) < 0.0)
			next.s = 1.0 - barelyMobilLeft(1.0 - next.s, props.s_oc);
		if ((1.0 - next.s - (1.0 - props.s_wc)) * (1.0 - iter.s - (1.0 - props.s_wc)) < 0.0)
			next.s = 1.0 - barelyMobilRight(1.0 - next.s, 1.0 - props.s_wc);
		// Water
		if ((next.s - props.s_wc) * (iter.s - props.s_wc) < 0.0)
			next.s = barelyMobilLeft(next.s, props.s_wc);
		if ((next.s - (1.0 - props.s_oc)) * (iter.s - (1.0 - props.s_oc)) < 0.0)
			next.s = barelyMobilRight(next.s, 1.0 - props.s_oc);
	};
	auto checkMaxResidual = [=, this](auto& next, auto& iter)
	{
		if (fabs(next.s - iter.s) > MAX_SAT_CHANGE)
			next.s = iter.s + sign(next.s - iter.s) * MAX_SAT_CHANGE;
	};

	for (size_t i = 0; i < size; i++)
	{
		auto& data = (*model)[i];
		checkCritPoints(data.u_next, data.u_iter, model->props_sk[0]);
		checkMaxResidual(data.u_next, data.u_iter);
	}
}
void Acid2dSolver::solveStep()
{
	int cellIdx, varIdx;
	err_newton = 1.0;
	averValue(averValPrev);
	std::fill(dAverVal.begin(), dAverVal.end(), 1.0);
	iterations = 0;

	auto continueIterations = [this]()
	{
		bool result = false;

		for (const auto& val : dAverVal)
			result += (val > CONV_VAR);
		
		return result * (err_newton > CONV_W2) * (iterations < MAX_ITER);
	};

	while (continueIterations())
	{
		copyIterLayer();

		computeJac();
		fill();
		solver.Assemble(ind_i, ind_j, a, elemNum, ind_rhs, rhs);
		solver.Solve(PRECOND::ILU_SIMPLE);
		copySolution(solver.getSolution());

		checkStability();
		err_newton = convergance(cellIdx, varIdx);

		averValue(averVal);
		for (int i = 0; i < var_size; i++)
			dAverVal[i] = fabs(averVal[i] - averValPrev[i]);
		averValPrev = averVal;

		model->snapshot_all(iterations + 1);
		iterations++;
	}

	cout << "Newton Iterations = " << iterations << endl;
}

void Acid2dSolver::computeJac()
{
	trace_on(0);

	for (size_t i = 0; i < size; i++)
	{
		model->x[i].m <<= model->u_next[var_size * i];
		model->x[i].p <<= model->u_next[var_size * i + 1];
		model->x[i].s <<= model->u_next[var_size * i + 2];
		model->x[i].xa <<= model->u_next[var_size * i + 3];
		model->x[i].xw <<= model->u_next[var_size * i + 4];
	}
	// Inner cells
	for (size_t i = 0; i < mesh->inner_cells; i++)
	{
		const auto& cell = mesh->cells[i];
		TapeVariable tmp = model->solveInner(cell);
		model->h[var_size * i] = tmp.m;
		model->h[var_size * i + 1] = tmp.p;
		model->h[var_size * i + 2] = tmp.s;
		model->h[var_size * i + 3] = tmp.xa;
		model->h[var_size * i + 4] = tmp.xw;
	}
	// Border cells
	for (size_t i = mesh->border_beg; i < model->cellsNum; i++)
	{
		const auto& cell = mesh->cells[i];
		TapeVariable tmp = model->solveBorder(cell);
		model->h[var_size * i] = tmp.m;
		model->h[var_size * i + 1] = tmp.p;
		model->h[var_size * i + 2] = tmp.s;
		model->h[var_size * i + 3] = tmp.xa;
		model->h[var_size * i + 4] = tmp.xw;
	}
/*	for (size_t i = 0; i < mesh->fracCells.size(); i++)
	{
		const auto& cell = *mesh->fracCells[i];
		model->h[var_size * cell.id] *= sqrt(cell.V);
		model->h[var_size * cell.id + 1] *= sqrt(cell.V);
		model->h[var_size * cell.id + 2] *= sqrt(cell.V);
		model->h[var_size * cell.id + 3] *= sqrt(cell.V);
		model->h[var_size * cell.id + 4] *= sqrt(cell.V);
	}*/
	// Well cell
	const int well_idx = mesh->well_idx;
	TapeVariable& cur = model->x[well_idx];
	model->h[well_idx * var_size + 1] = (cur.s - (1.0 - model->props_sk[0].s_oc)) / model->P_dim;
	model->h[well_idx * var_size + 2] += model->ht * model->props_w.getDensity(cur.p, cur.xa, cur.xw) * model->Q_sum / mesh->cells[well_idx].V;
	model->h[well_idx * var_size + 3] = (cur.xa - model->xa) / model->P_dim;
	model->h[well_idx * var_size + 4] = (cur.xw - (1.0 - model->xa)) / model->P_dim;

	for (size_t i = 0; i < size; i++)
	{
		for (size_t j = 0; j < var_size; j++)
			model->h[var_size * i + j] >>= y[var_size * i + j];
	}

	trace_off();
}
void Acid2dSolver::fill()
{
	sparse_jac(0, var_size * model->cellsNum, var_size * model->cellsNum, repeat,
		&model->u_next[0], &elemNum, (unsigned int**)(&ind_i), (unsigned int**)(&ind_j), &a, options);

	int counter = 0;
	for (const auto& cell : mesh->cells)
	{
		for (size_t i = 0; i < var_size; i++)
		{
			const int str_idx = var_size * cell.id + i;
			rhs[str_idx] = -y[str_idx];
		}
	}
}