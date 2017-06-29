#ifndef OIL2D_HPP_
#define OIL2D_HPP_

#include "src/models/Variables.hpp"
#include "src/models/AbstractModel.hpp"

#include "src/models/Oil2d/Properties.hpp"

namespace oil2d
{
	class Oil2d : public AbstractModel<var::containers::Var1phase, oil2d::Properties, var::BasicVariables, Oil2d>
	{
	protected:
		void setProps(Properties& props);
		void makeDimLess();
		void setPerforated();
		void setInitialState();

		std::vector<Skeleton_Props> props_sk;
		Oil_Props props_oil;
	public:
		Oil2d();
		~Oil2d();

		void setPeriod(const int period);
	};
};

#endif /* OIL2D_HPP_ */
