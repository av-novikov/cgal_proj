#ifndef ACID2D_HPP_
#define ACID2D_HPP_

#include <vector>

#include "src/models/Acid/Properties.hpp"
#include "src/models/Variables.hpp"
#include "src/models/AbstractModel.hpp"
#include "src/util/Interpolate.h"

namespace acid2d
{
	typedef CalciteReaction CurrentReaction;
	typedef CurrentReaction::REACTS REACTS;
	typedef var::containers::TapeAcidVar TapeVariable;

	class Acid2d : public AbstractModel<var::containers::AcidVar, Properties, var::BasicVariables, Acid2d>
	{
		template<typename> friend class VTKSnapshotter;
		friend class Acid2dSolver;
	protected: 
		CurrentReaction reac;
		Water_Props props_w;
		Oil_Props props_o;
		std::vector<Skeleton_Props> props_sk;

		TapeVariable* x;
		adouble* h;

		double xa;
		std::vector<double> xas;

		void setProps(const Properties& props);
		void makeDimLess();
		void setInitialState();

		// Service functions
		inline const size_t getUpwindIdx(size_t cur, size_t beta) const
		{
			if ((*this)[cur].u_next.p < (*this)[beta].u_next.p)
				return beta;
			else
				return cur;
		};
		inline adouble getReactionRate(const Cell& cell, const Skeleton_Props& props) const
		{
			const TapeVariable& var = x[cell.id];
			adouble isNotFrac = (cell.type == CellType::INNER || cell.type == CellType::BORDER) ? true : false;
			adouble tmp;
			condassign(tmp, isNotFrac, 
				var.s * props_w.getDensity(var.p, var.xa, var.xw) *	(var.xa - props.xa_eqbm) *
				reac.getReactionRate(props.m_init, var.m) / reac.comps[REACTS::ACID].mol_weight, 
				(adouble)0.0);
			return tmp;
		};
		const adouble getPerm(const Cell& cell) const
		{
			adouble isNotFrac = (cell.type == CellType::INNER || cell.type == CellType::BORDER) ? true : false;
			adouble tmp; 
			condassign(tmp, isNotFrac, props_sk[0].getPermCoseni_x(x[cell.id].m), (adouble)(props_sk[0].kx * 1000.0));
			return tmp;
		};
		double getPermValue(const Cell& cell) const
		{
			if (cell.type == CellType::INNER || cell.type == CellType::BORDER)
				return props_sk[0].getPermCoseni_x_value((*this)[cell.id].u_next.m);
			else 
				return props_sk[0].kx * 1000.0;
		};
		inline adouble getTrans(const Cell& cell, const size_t idx, const Cell& beta) const
		{
			adouble k1 = getPerm(cell);
			adouble k2 = getPerm(beta);
			const double dist1 = cell.dist[idx];
			const double dist2 = beta.getDistance(cell.id);
			return props_sk[0].height * cell.length[idx] * k1 * k2 / (k1 * dist2 + k2 * dist1);
		};

		TapeVariable solveInner(const Cell& cell);
		TapeVariable solveBorder(const Cell& cell);
	public:
		Acid2d();
		~Acid2d();

		void setPeriod(const int period);
		double getRate(const size_t cur);
		static const int var_size = VarContainer::size;
	};
};

#endif /* ACID2D_HPP_ */