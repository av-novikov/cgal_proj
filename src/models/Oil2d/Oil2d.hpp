#ifndef OIL2D_HPP_
#define OIL2D_HPP_

#include "src/models/Variables.hpp"
#include "src/models/AbstractModel.hpp"
#include "src/models/Oil2d/Properties.hpp"

namespace oil2d
{
	typedef var::containers::TapeVar1Phase TapeVariable;
	class Oil2d : public AbstractModel<var::containers::Var1phase, Properties, var::BasicVariables, Oil2d>
	{
		template<typename> friend class VTKSnapshotter;
		friend class Oil2dSolver;
	protected:
		void setProps(const Properties& props);
		void makeDimLess();
		void setInitialState();

		TapeVariable* x;
		adouble* h;

		std::vector<Skeleton_Props> props_sk;
		Oil_Props props_oil;

		const double getPerm(const Cell& cell) const
		{
			if (cell.type == CellType::INNER || cell.type == CellType::BORDER)
				return props_sk[0].kx;
			else
				return props_sk[0].kx * 1000.0;
		};
		double getTrans(const Cell& cell1, const int idx, const Cell& cell2) const
		{
			const double k1 = getPerm(cell1);
			const double k2 = getPerm(cell2);
			const double dist1 = cell1.dist[idx];
			const double dist2 = cell2.getDistance(cell1.id);
			return props_sk[0].height * cell1.length[idx] * k1 * k2 / (k1 * dist2 + k2 * dist1);
		};

		adouble solveInner(const Cell& cell);
		adouble solveBorder(const Cell& cell);
	public:
		Oil2d();
		~Oil2d();

		void setPeriod(const int period);
		static const int var_size = VarContainer::size;

		double getRate(const size_t cell_idx);
	};
};

#endif /* OIL2D_HPP_ */
