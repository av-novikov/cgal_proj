#ifndef OIL2DSOLVER_HPP_
#define OIL2DSOLVER_HPP_

#include "src/models/AbstractSolver.hpp"
#include "src/models/Oil2d/Oil2d.hpp"
#include "src/solvers/ParalutionInterface.h"
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
		ParSolver solver;

		void computeJac();
		void fill();
		void copySolution(const paralution::LocalVector<double>& sol);
	public:
		Oil2dSolver(Model* _model);
		~Oil2dSolver();

		void start();
	};
};

#endif /* OIL2DSOLVER_HPP_ */
