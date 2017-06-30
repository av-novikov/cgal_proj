#include "src/models/Oil2d/Oil2dSolver.hpp"
#include <iostream>

using namespace oil2d;
using namespace std;

Oil2dSolver::Oil2dSolver(Model* _model) : AbstractSolver<Model>(_model)
{
	plot_P.open("snaps/P.dat", ofstream::out);
	plot_Q.open("snaps/Q.dat", ofstream::out);
};
Oil2dSolver::~Oil2dSolver()
{
};
void Oil2dSolver::writeData()
{
	double p = 0.0, q = 0;

	plot_Q << cur_t * t_dim / 3600.0;

	map<int, double>::iterator it;
	for (it = model->Qcell.begin(); it != model->Qcell.end(); ++it)
	{
		p += model->u_next[it->first * Model::var_size] * model->P_dim;
		if (model->leftBoundIsRate)
			plot_Q << "\t" << it->second * model->Q_dim * 86400.0;
		else
		{
			plot_Q << "\t" << model->getRate(it->first) * model->Q_dim * 86400.0;
			q += model->getRate(it->first);
		}
	}
	plot_P << cur_t * t_dim / 3600.0 << "\t" << p / (double)(model->Qcell.size()) << endl;

	if (model->leftBoundIsRate)
		plot_Q << "\t" << model->Q_sum * model->Q_dim * 86400.0 << endl;
	else
		plot_Q << "\t" << q * model->Q_dim * 86400.0 << endl;
}
void Oil2dSolver::control()
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
void Oil2dSolver::solveStep()
{
	int cellIdx, varIdx;
	double err_newton = 1.0;
	double averPrev = averValue(0), aver, dAver = 1.0;

	iterations = 0;
	while (err_newton > 1.e-4 /*&& (dAverSat > 1.e-9 || dAverPres > 1.e-7)*/ && iterations < 20)
	{
		copyIterLayer();

		fill();
		//solver.Assemble(ind_i, ind_j, a, elemNum, ind_rhs, rhs);
		//solver.Solve();
		//copySolution(solver.getSolution());

		err_newton = convergance(cellIdx, varIdx);
		aver = averValue(0);		dAver = fabs(aver - averPrev);		averPrev = aver;
		iterations++;
	}

	cout << "Newton Iterations = " << iterations << endl;
}
void Oil2dSolver::copySolution(const paralution::LocalVector<double>& sol)
{
	for (int i = 0; i < model->cellsNum; i++)
	{
		auto& cell = (*model)[i];
		cell.u_next.p = sol[Model::var_size * i];
	}
}

void Oil2dSolver::fillIndices()
{
	int counter = 0;

	/*for (const auto& cell : mesh->cells)
	{
		getMatrixStencil(cell);
		for (const int idx : stencil_idx)
		{
			ind_i[counter] = Model::var_size * cell.id;			ind_j[counter++] = Model::var_size * idx;
		}
	}*/

	elemNum = counter;
	for (int i = 0; i < Model::var_size * model->cellsNum; i++)
		ind_rhs[i] = i;
}
void Oil2dSolver::fill()
{
}