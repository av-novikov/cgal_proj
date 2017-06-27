#include "src/models/AbstractSolver.hpp"
//#include "util/utils.h"

#include <iomanip>

using namespace std;

template <class modelType>
AbstractSolver<modelType>::AbstractSolver(modelType* _model) : model(_model), size(_model->getCellsNum()), Tt(model->period[model->period.size() - 1])
{
	NEWTON_STEP = 1.0;
	cur_t = cur_t_log = 0.0;
	curTimePeriod = 0;

	idx1 = int((model->perfIntervals[0].first + model->perfIntervals[0].second) / 2);
	//idx2 = idx1 + model->cellsNum_z + 1;

	t_dim = model->t_dim;
}
template <class modelType>
AbstractSolver<modelType>::~AbstractSolver()
{
}
template <class modelType>
void AbstractSolver<modelType>::start()
{
	int counter = 0;
	iterations = 8;

	model->setPeriod(curTimePeriod);
	while(cur_t < Tt)
	{
		control();
		if( model->isWriteSnaps )
			model->snapshot_all(counter++);
		doNextStep();
		copyTimeLayer();
		cout << "---------------------NEW TIME STEP---------------------" << endl;
		cout << setprecision(6);
		cout << "time = " << cur_t << endl;
	}
	if( model->isWriteSnaps )
		model->snapshot_all(counter++);
	writeData();
}
template <class modelType>
void AbstractSolver<modelType>::fill()
{
}
template <class modelType>
void AbstractSolver<modelType>::copyIterLayer()
{
	for (auto& cell : model->cells)
		cell.u_iter = cell.u_next;
}

template <class modelType>
void AbstractSolver<modelType>::revertIterLayer()
{
	for (int i = 0; i < model->cells.size(); i++)
		model->cells[i].u_next = model->cells[i].u_iter;
}
template <class modelType>
void AbstractSolver<modelType>::copyTimeLayer()
{
	for (auto& cell : model->cells)
		cell.u_prev = cell.u_iter = cell.u_next;
}

template <class modelType>
double AbstractSolver<modelType>::convergance(int& ind, int& varInd)
{
	double relErr = 0.0;
	double cur_relErr = 0.0;

	double var_next, var_iter;
		
	for(int i = 0; i < modelType::var_size; i++)
	{
		for(int j = 0; j < model->cells.size(); j++)
		{
			var_next = model->cells[j].u_next.values[i];	var_iter = model->cells[j].u_iter.values[i];
			if(fabs(var_next) > EQUALITY_TOLERANCE)	
			{
				cur_relErr = fabs( (var_next - var_iter) / var_next );
				if(cur_relErr > relErr)
				{
					relErr = cur_relErr;
					ind  = j;
					varInd = i;
				}
			}
		}
	}
	
	return relErr;
}
template <class modelType>
double AbstractSolver<modelType>::averValue(const int varInd)
{
	double tmp = 0.0;

	for(const auto& cell : model->cells)
	{
		tmp += cell.u_next.values[varInd] * cell.V;
	}
	
	return tmp / model->Volume;
}
template <class modelType>
void AbstractSolver<modelType>::averValue(std::array<double, modelType::var_size>& aver)
{
	std::fill(aver.begin(), aver.end(), 0.0);

	for (const auto& cell : model->cells)
		for (int i = 0; i < modelType::var_size; i++)
			aver[i] += cell.u_next.values[i] * cell.V;

	for(auto& val : aver)
		val /= model->Volume;
}

template <class modelType>
void AbstractSolver<modelType>::checkStability()
{
}
