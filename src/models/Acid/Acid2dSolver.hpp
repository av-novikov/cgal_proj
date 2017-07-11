#ifndef ACID2DSOLVER_HPP_
#define ACID2DSOLVER_HPP_

#include "src/models/AbstractSolver.hpp"
#include "src/models/Acid/Acid2d.hpp"
#include "src/solvers/ParalutionInterface.h"
#include <fstream>

namespace acid2d
{
	class Acid2dSolver : public AbstractSolver<Acid2d>
	{
	protected:
		void control();
		void solveStep();
		void writeData();

		std::array<double, var_size> averVal, averValPrev, dAverVal;
		double err_newton;

		std::ofstream S, P, qcells;
		ParSolver solver;

		void checkStability();
		void computeJac();
		void fill();
		void copySolution(const paralution::LocalVector<double>& sol);
	public:
		Acid2dSolver(acid2d::Acid2d* _model);
		~Acid2dSolver();

		void start();
	};
}

#endif /* ACID2DSOLVER_HPP_ */