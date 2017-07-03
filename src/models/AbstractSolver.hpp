#ifndef ABSTRACTSOLVER_HPP_
#define ABSTRACTSOLVER_HPP_

#include <algorithm>
#include <cmath>
#include <iostream>
#include <array>
#include <vector>

template <class modelType>
class AbstractSolver {
public:
	typedef modelType Model;
	typedef typename Model::Mesh Mesh;
	typedef typename Model::Cell Cell;
	typedef typename Model::TapeVariable TapeVariable;
protected:
	double** jac;
	double* y;

	Model* model;
	Mesh* mesh;
	const int size;
			
	int curTimePeriod;
	const double Tt;
	double cur_t, cur_t_log;

	int idx1, idx2;

	double t_dim;

	int iterations;
		
	void copyIterLayer();
	void revertIterLayer();
	void copyTimeLayer();
		
	double convergance(int& ind, int& varInd);
	double averValue(int varInd);
	void averValue(std::array<double, modelType::var_size>& aver);
		
	virtual void writeData() = 0;
	virtual void control() = 0;
	virtual void doNextStep();
	virtual void solveStep() = 0;
	double NEWTON_STEP;
	double CHOP_MULT;
	double MAX_SAT_CHANGE;
	double CONV_W2, CONV_VAR;
	int MAX_ITER;

	virtual void checkStability();

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
	};

	int* ind_i;
	int* ind_j;
	double* a;
	int* ind_rhs;
	double* rhs;
	// Number of non-zero elements in sparse matrix
	int elemNum;

public:
	AbstractSolver(modelType* _model);
	virtual ~AbstractSolver();
		
	virtual void fill();
	void fillIndices()
	{
		int counter = 0;

		for (const auto& cell : mesh->cells)
		{
			getMatrixStencil(cell);
			for (const int idx : stencil_idx)
			{
				ind_i[counter] = Model::var_size * cell.id;			ind_j[counter++] = Model::var_size * idx;
			}
			stencil_idx.clear();
		}

		elemNum = counter;
		for (int i = 0; i < Model::var_size * model->cellsNum; i++)
			ind_rhs[i] = i;
	}
	virtual void start();
};

#endif /* ABSTRACTSOLVER_HPP_ */
