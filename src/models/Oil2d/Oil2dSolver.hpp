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

		void fill();
		void fillIndices();
		void copySolution(const paralution::LocalVector<double>& sol);

		int* ind_i;
		int* ind_j;
		double* a;
		int* ind_rhs;
		double* rhs;
		// Number of non-zero elements in sparse matrix
		int elemNum;

		std::vector<int> stencil_idx;
		inline void getMatrixStencil(const Cell& cell)
		{
			if (cell.type == CellType::BORDER)
			{
				stencil_idx.resize(2);
				stencil_idx[0] = cell.id;
				stencil_idx[1] = cell.nebr[0];
			}
			else
			{
				stencil_idx.resize(4);
				stencil_idx[0] = cell.id;
				stencil_idx[1] = cell.nebr[0];
				stencil_idx[2] = cell.nebr[1];
				stencil_idx[3] = cell.nebr[2];
			}
		}
	public:
		Oil2dSolver(Model* _model);
		~Oil2dSolver();
	};
};

#endif /* OIL2DSOLVER_HPP_ */
