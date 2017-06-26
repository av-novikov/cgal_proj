#define _USE_MATH_DEFINES
#include <cmath>
#include "Model.hpp"

using namespace std;
using namespace point;
using namespace mesh;

int main(int argc, char* argv[])
{
	typedef Task::Body::Point Point;

	Task task;
	task.spatialStep = 0.5;
	Task::Body::Border body1border = { { 3, 3 },{ -3, 3 },{ -3, -3 },{ 3, -3 } };
	Task::Body::Border body2border = { { 3, 3 },{ 9, 3 },{ 9, -3 },{ 3, -3 } };
	task.bodies = {	Task::Body({ 0, body1border, {}})};

	Point pt1 = { -0.5, -1.0 };			Point pt2 = { -1.5, 0.5 };
	const int SIZE = 100;
	const double dx = (pt2[0] - pt1[0]) / (double)SIZE;
	const double dy = (pt2[1] - pt1[1]) / (double)SIZE;
	for (int i = 0; i < SIZE; i++)
	{
		Point p1 = {pt1[0] + (double)i * dx, pt1[1] + (double)i * dy};
		Point p2 = {pt1[0] + (double)(i+1) * dx, pt1[1] + (double)(i+1) * dy };
		task.bodies[0].constraint.push_back(make_pair(p1, p2));
	}
	
	FirstModel model;
	model.mesh = make_shared<FirstModel::Mesh>(task);
	model.mesh->snapshot(0);

	return 0;
}
