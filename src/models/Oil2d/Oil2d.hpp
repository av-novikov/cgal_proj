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
		typedef var::containers::TapeVar1Phase TapeVariables;
	protected:
		void setProps(Properties& props);
		void makeDimLess();
		void setPerforated();
		void setInitialState();

		std::vector<Skeleton_Props> props_sk;
		Oil_Props props_oil;

		void solveInner(const Cell& cell);
	public:
		Oil2d();
		~Oil2d();

		void setPeriod(const int period);
		static const int var_size = VarContainer::size;

		double getRate(const size_t cell_idx);
	};
};

#endif /* OIL2D_HPP_ */
