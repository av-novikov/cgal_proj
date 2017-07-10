#ifndef ACID2D_HPP_
#define ACID2D_HPP_

#include <vector>

#include "src/models/Acid/Properties.hpp"
#include "src/models/Variables.hpp"
#include "src/util/Interpolate.h"

namespace acid2d
{
	typedef AcidVar Variable;
	typedef TapeFirstAcid TapeVariable;
	typedef NewCylCell2D<Variable, Skeleton_Props> Cell;
	template <typename TVariable> using TCell = NewCylCell2D<TVariable, Skeleton_Props>;
	typedef CalciteReaction CurrentReaction;
	typedef CurrentReaction::REACTS REACTS;
	typedef Cell::Type Type;

	class Acid2d : public basic2d::Basic2d<Variable, Properties, Skeleton_Props, TCell, Acid2d>
	{
		template<typename> friend class Snapshotter;
		template<typename> friend class GRDECLSnapshotter;
		template<typename> friend class VTKSnapshotter;
		friend class Acid2dSolver;

	protected:
		CurrentReaction reac;
		Water_Props props_w;
		Oil_Props props_o;
		Gas_Props props_g;

		double xa;
		std::vector<double> xas;

		void setProps(Properties& props);
		void makeDimLess();
		void setInitialState();

		// Service functions
		inline adouble getAverage(adouble p1, const Cell& cell1, adouble p2, const Cell& cell2) const
		{
			double r1, r2;
			if (abs(cell1.num - cell2.num) == 1) 
			{
				r1 = cell1.hz;		r2 = cell2.hz;
			}
			else
			{
				r1 = cell1.hr;		r2 = cell2.hr;
			}

			return (p1 * (adouble)r2 + p2 * (adouble)r1) / (adouble)(r1 + r2);
		};
		inline adouble getReactionRate(TapeVariable& var, const Skeleton_Props& props) const
		{
			return var.sw * props_w.getDensity(var.p, var.xa, var.xw) *
					(var.xa - props.xa_eqbm) * 
					reac.getReactionRate(props.m_init, var.m) / reac.comps[REACTS::ACID].mol_weight;
		};
		inline adouble getTrans(const Cell& cell, adouble m_cell, const Cell& beta, adouble m_beta) const
		{
			adouble k1, k2, S;

			if (abs(cell.num - beta.num) == 1) {
				k1 = cell.props->getPermCoseni_z(m_cell);
				k2 = beta.props->getPermCoseni_z(m_beta);
				if (k1 == 0.0 && k2 == 0.0)
					return 0.0;
				S = 2.0 * M_PI * cell.r * cell.hr;
				return 2.0 * k1 * k2 * S / (k1 * beta.hz + k2 * cell.hz);
			}
			else {
				k1 = cell.props->getPermCoseni_r(m_cell);
				k2 = beta.props->getPermCoseni_r(m_beta);
				if (k1 == 0.0 && k2 == 0.0)
					return 0.0;
				S = 2.0 * M_PI * cell.hz * (cell.r + sign(beta.num - cell.num) * cell.hr / 2.0);
				return 2.0 * k1 * k2 * S / (k1 * beta.hr + k2 * cell.hr);
			}
		};
		inline double getWaterVelocity(Cell& cell, const int axis)
		{
			const Variable next = cell.u_next;

			switch (axis)
			{
			case R_AXIS:
				return -cell.props->getPermCoseni_r(next.m).value() * props_w.getKr(next.sw, next.so, cell.props).value() / props_w.getViscosity(next.p, next.xa, next.xw).value() * getNablaP(cell, NEXT, axis);
			case Z_AXIS:
				return -cell.props->getPermCoseni_z(next.m).value() * props_w.getKr(next.sw, next.so, cell.props).value() / props_w.getViscosity(next.p, next.xa, next.xw).value() * (getNablaP(cell, NEXT, axis) - props_w.getDensity(next.p, next.xa, next.xw).value() * grav);
			}
		};
		inline double getGasVelocity(Cell& cell, const int axis)
		{
			const Variable next = cell.u_next;

			switch (axis)
			{
			case R_AXIS:
				return -cell.props->getPermCoseni_r(next.m).value() * props_g.getKr(next.sw, next.so, cell.props).value() / props_g.getViscosity(next.p).value() * getNablaP(cell, NEXT, axis);
			case Z_AXIS:
				return -cell.props->getPermCoseni_z(next.m).value() * props_g.getKr(next.sw, next.so, cell.props).value() / props_g.getViscosity(next.p).value() * (getNablaP(cell, NEXT, axis) - props_g.getDensity(next.p).value() * grav);
			}
		};
		inline double getOilVelocity(Cell& cell, const int axis)
		{
			const Variable next = cell.u_next;

			switch (axis)
			{
			case R_AXIS:
				return -cell.props->getPermCoseni_r(next.m).value() * props_o.getKr(next.sw, next.so, cell.props).value() / props_o.getViscosity(next.p).value() * getNablaP(cell, NEXT, axis);
			case Z_AXIS:
				return -cell.props->getPermCoseni_z(next.m).value() * props_o.getKr(next.sw, next.so, cell.props).value() / props_o.getViscosity(next.p).value() * (getNablaP(cell, NEXT, axis) - props_o.getDensity(next.p, next.p_bub, next.SATUR).value() * grav);
			}
		};

		void solve_eqMiddle(const Cell& cell);
		void solve_eqLeft(const Cell& cell);
		void solve_eqRight(const Cell& cell);
		void solve_eqVertical(const Cell& cell);
		void setVariables(const Cell& cell);

	public:
		Acid2d();
		~Acid2d();

		double* x;
		double* y;
		double** jac;

		void setPeriod(int period);
		double getRate(int cur);
		static const int var_size = Variable::size - 1;
	};
};

#endif /* ACID2D_HPP_ */