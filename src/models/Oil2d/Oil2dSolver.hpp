#ifndef OIL2DSOLVER_HPP_
#define OIL2DSOLVER_HPP_

#include "src/models/AbstractSolver.hpp"
#include "src/models/Oil2d/Oil2d.hpp"
#include <fstream>

namespace oil2d
{
	class Oil2dSolver : public AbstractSolver<Oil2d>
	{
	protected:
		void control();
		void solveStep();
		void writeData();

		std::ofstream plot_P, plot_Q;
	public:
		Oil2dSolver(Model* _model);
		~Oil2dSolver();
	};
};

#endif /* OIL2DSOLVER_HPP_ */
