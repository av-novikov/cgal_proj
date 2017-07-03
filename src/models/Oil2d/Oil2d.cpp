#include "src/models/Oil2d/Oil2d.hpp"
#include "src/util/utils.h"

#include "adolc/drivers/drivers.h"
#include "adolc/adolc.h"

using namespace oil2d;
using namespace std;

Oil2d::Oil2d()
{
}
Oil2d::~Oil2d()
{

	delete[] x, h;
}
void Oil2d::setProps(const Properties& props)
{
	r_w = props.r_w;
	r_e = props.r_e;
	perfIntervals = props.perfIntervals;

	leftBoundIsRate = props.leftBoundIsRate;
	rightBoundIsPres = props.rightBoundIsPres;

	skeletonsNum = props.props_sk.size();
	props_sk = props.props_sk;
	for (int j = 0; j < skeletonsNum; j++)
	{
		props_sk[j].kx = MilliDarcyToM2(props_sk[j].kx);
		props_sk[j].ky = MilliDarcyToM2(props_sk[j].ky);
	}

	periodsNum = props.timePeriods.size();
	for (int i = 0; i < periodsNum; i++)
	{
		period.push_back(props.timePeriods[i]);

		if (leftBoundIsRate)
			rate.push_back(props.rates[i] / 86400.0);
		else
			pwf.push_back(props.pwf[i]);
	}

	ht = props.ht;
	ht_min = props.ht_min;
	ht_max = props.ht_max;

	props_oil = props.props_oil;
	props_oil.visc = cPToPaSec(props_oil.visc);

	alpha = props.alpha;

	makeDimLess();
}
void Oil2d::makeDimLess()
{
	R_dim = r_w;
	t_dim = 3600.0;
	P_dim = props_sk[0].p_init;

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

	props_oil.visc /= (P_dim * t_dim);
	props_oil.dens_stc /= (P_dim * t_dim * t_dim / R_dim / R_dim);
	props_oil.beta /= (1.0 / P_dim);

	alpha /= t_dim;
}
void Oil2d::setPerforated()
{
	// Well locating & preparing
	Qcell[0] = 0.0;

	x = new TapeVariable[cellsNum];
	h = new adouble[var_size * cellsNum];
}
void Oil2d::setInitialState()
{
	const auto& props = props_sk[0];
	for (size_t i = 0; i < cellsNum; i++)
	{
		auto cell = (*this)[i];
		cell.u_prev.p = cell.u_iter.p = cell.u_next.p = props.p_init;
	}
}
void Oil2d::setPeriod(const int period)
{
	if (leftBoundIsRate)
	{
		Q_sum = rate[period];

		/*if (period == 0 || rate[period - 1] < EQUALITY_TOLERANCE) {
			map<int, double>::iterator it;
			for (it = Qcell.begin(); it != Qcell.end(); ++it)
				it->second = Q_sum * cells[it->first].hz / height_perf;
		}
		else {
			map<int, double>::iterator it;
			for (it = Qcell.begin(); it != Qcell.end(); ++it)
				it->second = it->second * Q_sum / rate[period - 1];
		}*/
	}
	else
	{
		Pwf = pwf[period];
		Q_sum = 0.0;
	}
}
double Oil2d::getRate(const size_t cell_idx)
{
	return 0.0;
}

adouble Oil2d::solveInner(const Cell& cell)
{
	const auto& cur = x[cell.id];
	const auto& prev = (*this)[cell.id].u_prev;

	adouble H = props_sk[0].getPoro(cur.p) * props_sk[0].getDensity(cur.p) - props_sk[0].getPoro(prev.p) * props_sk[0].getDensity(prev.p);
	for (const int nebr_idx : cell.nebr)
	{
		const auto& beta = mesh->cells[nebr_idx];
		const auto& nebr = x[nebr_idx];

		H += ht / cell.V / props_oil.getViscosity(nebr.p);
	}
	return H;
}
adouble Oil2d::solveBorder(const Cell& cell)
{
	return 0;
}