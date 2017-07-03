#ifndef OIL2D_HPP_
#define OIL2D_HPP_

#include "src/models/Variables.hpp"
#include "src/models/AbstractModel.hpp"

#include "src/models/Oil2d/Properties.hpp"

namespace oil2d
{
	class Oil2d : public AbstractModel<var::containers::Var1phase, oil2d::Properties, var::BasicVariables, Oil2d>
	{
		friend class Oil2dSolver;
	public:
		typedef var::containers::TapeVar1Phase TapeVariable;
	protected:
		void setProps(const Properties& props);
		void makeDimLess();
		void setPerforated();
		void setInitialState();

		TapeVariable* x;
		adouble* h;

		std::vector<Skeleton_Props> props_sk;
		Oil_Props props_oil;

		double getPerm(const Cell& cell)
		{
			if (cell.type == CellType::INNER)
				return props_sk[0].kx;
			else if (cell.type == CellType::FRAC)
				return props_sk[0].kx * 100.0;
		};
		/*double getTrans(const Cell& cell1, const Cell& cell2)
		{
			const double k1 = getPerm(cell1);
			const double k2 = getPerm(cell2);
			const double S = props_sk[0].height * props
		};*/

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
