#include "src/models/Acid/Acid2d.hpp"

#include <cassert>

#include "adolc/drivers/drivers.h"
#include "adolc/adolc.h"

using namespace acid2d;

double acid2d::Component::R = 8.3144598;
double acid2d::Component::p_std = 101325.0;

Acid2d::Acid2d()
{
}
Acid2d::~Acid2d()
{
	delete[] x, h;
}
void Acid2d::setProps(Properties& props)
{
	R_dim = props.R_dim;
	r_w = props.r_w;
	r_e = props.r_e;

	leftBoundIsRate = props.leftBoundIsRate;
	rightBoundIsPres = props.rightBoundIsPres;

	perfIntervals = props.perfIntervals;

	props_sk = props.props_sk;
	skeletonsNum = props.props_sk.size();
	for (int j = 0; j < skeletonsNum; j++)
	{
		props_sk[j].kx = MilliDarcyToM2(props_sk[j].kx);
		props_sk[j].ky = MilliDarcyToM2(props_sk[j].ky);
	}

	periodsNum = props.timePeriods.size();
	for (int i = 0; i < periodsNum; i++)
	{
		period.push_back(props.timePeriods[i]);
		xas.push_back(props.xa[i]);
		if (leftBoundIsRate)
			rate.push_back(props.rates[i] / 86400.0);
		else
			pwf.push_back(props.pwf[i]);
	}

	ht = props.ht;
	ht_min = props.ht_min;
	ht_max = props.ht_max;

	props_w = props.props_w;
	props_w.visc = cPToPaSec(props_w.visc);

	props_o = props.props_o;
	props_o.visc = cPToPaSec(props_o.visc);

	for (auto& comp : reac.comps)
		comp.mol_weight = gramToKg(comp.mol_weight);

	makeDimLess();

	// Data sets
	//props_o.b = setDataset(props.B_oil, P_dim / BAR_TO_PA, 1.0);
	//props_o.Rs = setDataset(props.Rs, P_dim / BAR_TO_PA, 1.0);
	//props_g.rho = setDataset(props.rho_co2, P_dim / BAR_TO_PA, (P_dim * t_dim * t_dim / R_dim / R_dim));
}
void Acid2d::makeDimLess()
{
	T_dim = Component::T;
	t_dim = 3600.0;
	P_dim = props_sk[0].p_init;

	Component::p_std /= P_dim;
	Component::R /= (P_dim * R_dim * R_dim * R_dim / T_dim);
	Component::T /= T_dim;

	// Temporal properties
	ht /= t_dim;
	ht_min /= t_dim;
	ht_max /= t_dim;

	// Grid properties
	r_w /= R_dim;
	r_e /= R_dim;

	for (auto& props : props_sk)
	{
		props.kx /= (R_dim * R_dim);
		props.ky /= (R_dim * R_dim);
		props.dens_stc /= (P_dim * t_dim * t_dim / R_dim / R_dim);
		props.beta /= (1.0 / P_dim);
		props.height /= R_dim;
		props.p_init /= P_dim;
		props.p_out /= P_dim;
	}

	Q_dim = R_dim * R_dim * R_dim / t_dim;
	for (int i = 0; i < periodsNum; i++)
	{
		period[i] /= t_dim;
		if (leftBoundIsRate)
			rate[i] /= Q_dim;
		else
			pwf[i] /= P_dim;
	}

	props_w.visc /= (P_dim * t_dim);
	props_w.dens_stc /= (P_dim * t_dim * t_dim / R_dim / R_dim);
	props_w.beta /= (1.0 / P_dim);
	props_o.visc /= (P_dim * t_dim);
	props_o.dens_stc /= (P_dim * t_dim * t_dim / R_dim / R_dim);
	props_o.gas_dens_stc /= (P_dim * t_dim * t_dim / R_dim / R_dim);
	props_o.beta /= (1.0 / P_dim);

	reac.activation_energy /= (P_dim * R_dim * R_dim * R_dim);
	reac.surf_init /= (1.0 / R_dim);
	reac.reaction_const /= (R_dim / t_dim);
	for (auto& comp : reac.comps)
	{
		comp.rho_stc /= (P_dim * t_dim * t_dim / R_dim / R_dim);
		comp.mol_weight /= (P_dim * t_dim * t_dim * R_dim);
	}
}
void Acid2d::setPeriod(int period)
{
	if (leftBoundIsRate)
	{
		Q_sum = rate[period];
	}
	else
	{
		Pwf = pwf[period];
		Q_sum = 0.0;
	}

	xa = xas[period];
}
void Acid2d::setInitialState()
{
	const auto& props = props_sk[0];
	for (size_t i = 0; i < cellsNum; i++)
	{
		const auto& cell = mesh->cells[i];
		auto data = (*this)[i];

		if (cell.type == CellType::FRAC)
			data.u_prev.p = data.u_iter.p = data.u_next.p = props.p_init * 0.99;
		else
			data.u_prev.p = data.u_iter.p = data.u_next.p = props.p_init;
		data.u_prev.m = data.u_iter.m = data.u_next.m = props.m_init;
		data.u_prev.s = data.u_iter.s = data.u_next.s = props.s_init;
		data.u_prev.xa = data.u_iter.xa = data.u_next.xa = props.xa_init;
		data.u_prev.xw = data.u_iter.xw = data.u_next.xw = props.xw_init;
	}

	x = new TapeVariable[cellsNum];
	h = new adouble[var_size * cellsNum];
}
double Acid2d::getRate(int cur)
{
	return 0.0;
};

TapeVariable Acid2d::solveInner(const Cell& cell)
{
	const auto& cur = x[cell.id];
	const auto& prev = (*this)[cell.id].u_prev;
	const auto& props = props_sk[0];

	adouble rate = getReactionRate(cur, props);
	TapeVariable res;
	res.m = (1.0 - cur.m) * props.getDensity(cur.p) - (1.0 - prev.m) * props.getDensity(prev.p) -
		ht * reac.indices[REACTS::CALCITE] * reac.comps[REACTS::CALCITE].mol_weight * rate;
	res.p = cur.m * cur.s * props_w.getDensity(cur.p, cur.xa, cur.xw) -
		prev.m * prev.s * props_w.getDensity(prev.p, prev.xa, prev.xw) -
		ht * (reac.indices[REACTS::ACID] * reac.comps[REACTS::ACID].mol_weight +
			reac.indices[REACTS::WATER] * reac.comps[REACTS::WATER].mol_weight +
			reac.indices[REACTS::SALT] * reac.comps[REACTS::SALT].mol_weight) * rate;
	res.s = cur.m * cur.s * props_w.getDensity(cur.p, cur.xa, cur.xw) * cur.xa -
		prev.m * prev.s * props_w.getDensity(prev.p, prev.xa, prev.xw) * prev.xa -
		ht * reac.indices[REACTS::ACID] * reac.comps[REACTS::ACID].mol_weight * rate;
	res.xa = cur.m * cur.s * props_w.getDensity(cur.p, cur.xa, cur.xw) * cur.xw -
		prev.m * prev.s * props_w.getDensity(prev.p, prev.xa, prev.xw) * prev.xw -
		ht * reac.indices[REACTS::WATER] * reac.comps[REACTS::WATER].mol_weight * rate;
	H[4] = cur.m * cur.s * props_o.dens_stc / props_o.getB(cur.p) -
		prev.m * prev.s * props_o.dens_stc / props_o.getB(prev.p);

	for (int i = 0; i < 3; i++)
	{
		const int nebr_idx = cell.nebr[i];
		const auto& beta = mesh->cells[nebr_idx];
		const auto& nebr = x[nebr_idx];
		const int upwd_idx = getUpwindIdx(cell.id, beta.id);
		TapeVariable& upwd = x[upwd_idx];

		adouble dens_w = linearAppr(props_w.getDensity(cur.p, cur.xa, cur.xw), cell.dist[i], 
									props_w.getDensity(nebr.p, nebr.xa, nebr.xw), beta.getDistance(cell.id));
		adouble dens_o = linearAppr(props_o.getDensity(cur.p), cell.dist[i],
									props_o.getDensity(nebr.p), beta.getDistance(cell.id));
		adouble buf_w = ht / cell.V * getTrans(cell, i, beta) * (cur.p - nebr.p) *
			dens_w * props_w.getKr(upwd.s, props) / props_w.getViscosity(upwd.p, upwd.xa, upwd.xw);
		adouble buf_o = ht / cell.V * getTrans(cell, i, beta) * (cur.p - nebr.p) *
			dens_o * props_o.getKr(upwd.s, props) / props_o.getViscosity(upwd.p);

		h[1] += buf_w;
		h[2] += buf_w * upwd.xa;
		h[3] += buf_w * upwd.xw;
		h[4] += buf_o;
	}
	return res;
}

void Acid2d::solve_eqMiddle(const Cell& cell)
{
	const Skeleton_Props& props = *cell.props;

	trace_on(mid);
	adouble h[var_size], tmp;
	TapeVariable var[stencil];
	const Variable& prev = cell.u_prev;

	for (int i = 0; i < stencil; i++)
	{
		var[i].m <<= x[i * Variable::size];
		var[i].p <<= x[i * Variable::size + 1];
		var[i].sw <<= x[i * Variable::size + 2];
		var[i].xa <<= x[i * Variable::size + 3];
		var[i].xw <<= x[i * Variable::size + 4];
		var[i].so <<= x[i * Variable::size + 5];
		var[i].p_bub <<= x[i * Variable::size + 6];
	}

	TapeVariable& next = var[0];
	adouble satur = cell.u_next.SATUR;
	adouble rate = getReactionRate(next, props);
	double rate0 = getReactionRate(next, props).value();

	h[0] = (1.0 - next.m) * props.getDensity(next.p) -
			(1.0 - prev.m) * props.getDensity(prev.p) - 
			ht * reac.indices[REACTS::CALCITE] * reac.comps[REACTS::CALCITE].mol_weight * rate;
	/*h[1] = next.m * (1.0 - next.sw - next.so) * props_g.getDensity(next.p) -
		prev.m * (1.0 - prev.sw - prev.so) * props_g.getDensity(prev.p) -
		ht * reac.indices[REACTS::CO2] * reac.comps[REACTS::CO2].mol_weight * rate;*/
	condassign(h[1], satur,
				next.m * ((1.0 - next.sw - next.so) * props_g.getDensity(next.p) +
			next.so * props_o.gas_dens_stc * props_o.getRs(next.p, next.p_bub, satur) / props_o.getB(next.p, next.p_bub, satur)) -
				prev.m * ((1.0 - prev.sw - prev.so) * props_g.getDensity(prev.p) + 
			prev.so * props_o.gas_dens_stc * props_o.getRs(prev.p, prev.p_bub, prev.SATUR) / props_o.getB(prev.p, prev.p_bub, prev.SATUR)) -
				ht * reac.indices[REACTS::CO2] * reac.comps[REACTS::CO2].mol_weight * rate,
				next.m * next.so * props_o.gas_dens_stc * props_o.getRs(next.p, next.p_bub, satur) / props_o.getB(next.p, next.p_bub, satur) -
			prev.m * prev.so * props_o.gas_dens_stc * props_o.getRs(prev.p, prev.p_bub, prev.SATUR) / props_o.getB(prev.p, prev.p_bub, prev.SATUR) - 
				ht * reac.indices[REACTS::CO2] * reac.comps[REACTS::CO2].mol_weight * rate);
	h[2] = next.m * next.sw * props_w.getDensity(next.p, next.xa, next.xw) -
		prev.m * prev.sw * props_w.getDensity(prev.p, prev.xa, prev.xw) -
		ht * (reac.indices[REACTS::ACID] * reac.comps[REACTS::ACID].mol_weight +
			reac.indices[REACTS::WATER] * reac.comps[REACTS::WATER].mol_weight +
			reac.indices[REACTS::SALT] * reac.comps[REACTS::SALT].mol_weight) * rate;
	h[3] = next.m * next.sw * props_w.getDensity(next.p, next.xa, next.xw) * next.xa -
		prev.m * prev.sw * props_w.getDensity(prev.p, prev.xa, prev.xw) * prev.xa -
		ht * reac.indices[REACTS::ACID] * reac.comps[REACTS::ACID].mol_weight * rate;
	h[4] = next.m * next.sw * props_w.getDensity(next.p, next.xa, next.xw) * next.xw -
		prev.m * prev.sw * props_w.getDensity(prev.p, prev.xa, prev.xw) * prev.xw -
		ht * reac.indices[REACTS::WATER] * reac.comps[REACTS::WATER].mol_weight * rate;
	h[5] = next.m * next.so * props_o.dens_stc / props_o.getB(next.p, next.p_bub, satur) -
			prev.m * prev.so * props_o.dens_stc / props_o.getB(prev.p, prev.p_bub, prev.SATUR);

	int neighbor[4];
	getNeighborIdx(cell.num, neighbor);
	for (int i = 0; i < 4; i++)
	{
		const Cell& beta = cells[neighbor[i]];
		const int upwd_idx = (getUpwindIdx(cell.num, neighbor[i]) == cell.num) ? 0 : i + 1;
		const TapeVariable& nebr = var[i + 1];
		TapeVariable& upwd = var[upwd_idx];
		
		adouble dens_w = getAverage(props_w.getDensity(next.p, next.xa, next.xw), cell, props_w.getDensity(nebr.p, nebr.xa, nebr.xw), beta);
		adouble dens_o = getAverage(props_o.dens_stc / props_o.getB(next.p, next.p_bub, satur), cell, 
									props_o.dens_stc / props_o.getB(nebr.p, nebr.p_bub, beta.u_next.SATUR), beta);
		adouble dens_g = getAverage(props_g.getDensity(next.p), cell, props_g.getDensity(nebr.p), beta);
		adouble dens_go = props_o.gas_dens_stc *
				getAverage(props_o.getRs(next.p, next.p_bub, satur) / props_o.getB(next.p, next.p_bub, satur), cell,
							props_o.getRs(nebr.p, nebr.p_bub, beta.u_next.SATUR) / props_o.getB(nebr.p, nebr.p_bub, beta.u_next.SATUR), beta);
		
		adouble buf_g = ht / cell.V * getTrans(cell, next.m, beta, nebr.m) * ((next.p - nebr.p) - dens_g * grav * (cell.z - beta.z)) *
						dens_g * props_g.getKr(upwd.sw, upwd.so, cells[upwd_idx].props) / props_g.getViscosity(upwd.p);
		adouble buf_w = ht / cell.V * getTrans(cell, next.m, beta, nebr.m) * ((next.p - nebr.p) - dens_w * grav * (cell.z - beta.z)) *
						dens_w * props_w.getKr(upwd.sw, upwd.so, cells[upwd_idx].props) / props_w.getViscosity(upwd.p, upwd.xa, upwd.xw);
		adouble buf_go = ht / cell.V * getTrans(cell, next.m, beta, nebr.m) * ((next.p - nebr.p) - dens_o * grav * (cell.z - beta.z)) *
						dens_go * props_o.getKr(upwd.sw, upwd.so, cells[upwd_idx].props) / props_o.getViscosity(upwd.p);
		adouble buf_o = ht / cell.V * getTrans(cell, next.m, beta, nebr.m) * ((next.p - nebr.p) - dens_o * grav * (cell.z - beta.z)) *
			dens_o * props_o.getKr(upwd.sw, upwd.so, cells[upwd_idx].props) / props_o.getViscosity(upwd.p);

		condassign(tmp, satur, buf_g + buf_go, buf_go);
		h[1] += tmp;
		h[2] += buf_w;
		h[3] += buf_w * upwd.xa;
		h[4] += buf_w * upwd.xw;
		h[5] += buf_o;
	}

	//for (int i = 0; i < var_size; i++)
	//	h[i] /= P_dim;

	for (int i = 0; i < var_size; i++)
		h[i] >>= y[i];

	trace_off();
}
void Acid2d::solve_eqLeft(const Cell& cell)
{
	const Cell& beta1 = cells[cell.num + cellsNum_z + 2];
	const Cell& beta2 = cells[cell.num + 2 * cellsNum_z + 4];

	trace_on(left);
	adouble h[var_size];
	TapeVariable var[Lstencil];
	for (int i = 0; i < Lstencil; i++)
	{
		var[i].m <<= x[i * Variable::size];
		var[i].p <<= x[i * Variable::size + 1];
		var[i].sw <<= x[i * Variable::size + 2];
		var[i].xa <<= x[i * Variable::size + 3];
		var[i].xw <<= x[i * Variable::size + 4];
		var[i].so <<= x[i * Variable::size + 5];
		var[i].p_bub <<= x[i * Variable::size + 6];
	}

	const adouble leftIsRate = leftBoundIsRate;
	TapeVariable& next = var[0];
	TapeVariable& nebr1 = var[1];
	TapeVariable& nebr2 = var[2];
	const Variable& prev = cell.u_prev;
	const Skeleton_Props& props = *cell.props;

	adouble isPerforated, leftPresPerforated, leftPresNonPerforated;
	auto it = Qcell.find(cell.num);
	double rate;
	if (it != Qcell.end())
	{
		leftPresPerforated = !leftBoundIsRate;
		leftPresNonPerforated = false;
		isPerforated = true;
		rate = Qcell[cell.num];
	}
	else
	{
		isPerforated = false;
		leftPresPerforated = false;
		leftPresNonPerforated = !leftBoundIsRate;
		rate = 0.0;
	}

	adouble nonsatur = !cell.u_next.SATUR;
	condassign(h[0], isPerforated, (1.0 - next.m) * props.getDensity(next.p) - (1.0 - prev.m) * props.getDensity(prev.p) -
		ht * reac.indices[REACTS::CALCITE] * reac.comps[REACTS::CALCITE].mol_weight * getReactionRate(next, props),
									next.m - nebr1.m);
	
	condassign(h[1], leftIsRate, props_w.getDensity(next.p, next.xa, next.xw) * getTrans(cell, next.m, beta1, nebr1.m) /
		props_w.getViscosity(next.p, next.xa, next.xw) * (nebr1.p - next.p) +
		props_w.getDensity(Component::p_std, next.xa, next.xw) * rate,
									next.p - Pwf);
	condassign(h[1], leftPresPerforated, next.p - Pwf);
	condassign(h[1], leftPresNonPerforated, next.p - nebr1.p);

	condassign(h[2], isPerforated, next.sw - (1.0 - props.s_oc), 
									next.sw - nebr1.sw);
	condassign(h[3], isPerforated, next.xa - xa,
									next.xa - nebr1.xa);
	condassign(h[4], isPerforated, next.xw - (1.0 - xa), 
									next.xw - nebr1.xw);
	condassign(h[5], isPerforated, next.so - props.s_oc,
		next.so - nebr1.so);
	//condassign(h[5], nonsatur, next.p_bub - nebr.p_bub);
	condassign(h[5], nonsatur, (next.p_bub - nebr1.p_bub) / (adouble)(cell.r - beta1.r) -
								(nebr1.p_bub - nebr2.p_bub) / (adouble)(beta1.r - beta2.r));

	//for (int i = 0; i < var_size; i++)
	//	h[i] /= P_dim;

	for (int i = 0; i < var_size; i++)
		h[i] >>= y[i];

	trace_off();
}
void Acid2d::solve_eqRight(const Cell& cell)
{
	trace_on(right);
	adouble h[var_size];
	TapeVariable var[Rstencil];
	adouble rightIsPres = rightBoundIsPres;
	for (int i = 0; i < Rstencil; i++)
	{
		var[i].m <<= x[i * Variable::size];
		var[i].p <<= x[i * Variable::size + 1];
		var[i].sw <<= x[i * Variable::size + 2];
		var[i].xa <<= x[i * Variable::size + 3];
		var[i].xw <<= x[i * Variable::size + 4];
		var[i].so <<= x[i * Variable::size + 5];
		var[i].p_bub <<= x[i * Variable::size + 6];
	}

	TapeVariable& next = var[0];
	TapeVariable& nebr = var[1];
	const Variable& prev = cell.u_prev;
	const Skeleton_Props& props = *cell.props;
	adouble satur = cell.u_next.SATUR;

	h[0] = next.m - nebr.m;
	condassign(h[1], rightIsPres, next.p - (adouble)(cell.props->p_out), next.p - (adouble)(nebr.p));
	h[2] = next.sw - nebr.sw;
	h[3] = next.xa - nebr.xa;
	h[4] = next.xw - nebr.xw;
	condassign(h[5], satur, next.so - nebr.so, next.p_bub - nebr.p_bub);
	
	//for (int i = 0; i < var_size; i++)
	//	h[i] /= P_dim;

	for (int i = 0; i < var_size; i++)
		h[i] >>= y[i];

	trace_off();
}
void Acid2d::solve_eqVertical(const Cell& cell)
{
	trace_on(vertical);
	adouble h[var_size];
	TapeVariable var[Vstencil];

	for (int i = 0; i < Vstencil; i++)
	{
		var[i].m <<= x[i * Variable::size];
		var[i].p <<= x[i * Variable::size + 1];
		var[i].sw <<= x[i * Variable::size + 2];
		var[i].xa <<= x[i * Variable::size + 3];
		var[i].xw <<= x[i * Variable::size + 4];
		var[i].so <<= x[i * Variable::size + 5];
		var[i].p_bub <<= x[i * Variable::size + 6];
	}

	const TapeVariable& next = var[0];
	const TapeVariable& nebr = var[1];
	adouble satur = cell.u_next.SATUR;

	h[0] = next.m - nebr.m;
	h[1] = next.p - nebr.p;
	h[2] = next.sw - nebr.sw;
	h[3] = next.xa - nebr.xa;
	h[4] = next.xw - nebr.xw;
	condassign(h[5], satur, next.so - nebr.so, next.p_bub - nebr.p_bub);

	//for (int i = 0; i < var_size; i++)
	//	h[i] /= P_dim;

	for (int i = 0; i < var_size; i++)
		h[i] >>= y[i];

	trace_off();
}
void Acid2d::setVariables(const Cell& cell)
{
	if (cell.type == Type::WELL_LAT)
	{
		// Left
		const Variable& next = cell.u_next;
		const Variable& nebr1 = cells[cell.num + cellsNum_z + 2].u_next;
		const Variable& nebr2 = cells[cell.num + 2 * cellsNum_z + 4].u_next;
		for (int i = 0; i < Variable::size; i++)
		{
			x[i] = next.values[i];
			x[Variable::size + i] = nebr1.values[i];
			x[2 * Variable::size + i] = nebr2.values[i];
		}

		solve_eqLeft(cell);
		jacobian(left, var_size, Variable::size * Lstencil, x, jac);
	}
	else if (cell.type == Type::RIGHT)
	{
		// Right
		const Variable& next = cells[cell.num].u_next;
		const Variable& nebr1 = cells[cell.num - cellsNum_z - 2].u_next;
		const Variable& nebr2 = cells[cell.num - 2 * cellsNum_z - 4].u_next;
		for (int i = 0; i < Variable::size; i++)
		{
			x[i] = next.values[i];
			x[Variable::size + i] = nebr1.values[i];
			x[2 * Variable::size + i] = nebr2.values[i];
		}

		solve_eqRight(cell);
		jacobian(right, var_size, Variable::size * Rstencil, x, jac);
	}
	else if (cell.type == Type::TOP)
	{
		// Top
		const Variable& next = cells[cell.num].u_next;
		const Variable& nebr = cells[cell.num + 1].u_next;
		for (int i = 0; i < Variable::size; i++)
		{
			x[i] = next.values[i];
			x[Variable::size + i] = nebr.values[i];
		}

		solve_eqVertical(cell);
		jacobian(vertical, var_size, Variable::size * Vstencil, x, jac);
	}
	else if (cell.type == Type::BOTTOM)
	{
		// Bottom
		const Variable& next = cells[cell.num].u_next;
		const Variable& nebr = cells[cell.num - 1].u_next;
		for (int i = 0; i < Variable::size; i++)
		{
			x[i] = next.values[i];
			x[Variable::size + i] = nebr.values[i];
		}

		solve_eqVertical(cell);
		jacobian(vertical, var_size, Variable::size * Vstencil, x, jac);
	}
	else if (cell.type == Type::MIDDLE)
	{
		// Middle
		const Variable& next = cells[cell.num].u_next;
		int neighbor[stencil - 1];
		getNeighborIdx(cell.num, neighbor);
		for (int i = 0; i < Variable::size; i++)
		{
			x[i] = next.values[i];

			for (int j = 0; j < stencil - 1; j++)
			{
				const Variable& nebr = cells[neighbor[j]].u_next;
				x[(j + 1) * Variable::size + i] = nebr.values[i];
			}
		}

		solve_eqMiddle(cell);
		jacobian(mid, var_size, Variable::size * stencil, x, jac);
	}
}