#define _USE_MATH_DEFINES
#include <cmath>
#include "src/models/Oil2d/Oil2d.hpp"

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

	const double w = 0.01;
	Point pt1 = { -0.5, -1.0 };			Point pt2 = { -1.5, 0.5 };
	Point pt3 = pt2;		pt3[0] += w;
	Point pt4 = pt1;		pt4[0] += w;
	
	const int SIZE = 1;
	double dx = (pt2[0] - pt1[0]) / (double)SIZE;
	double dy = (pt2[1] - pt1[1]) / (double)SIZE;
	Point p1, p2;
	for (int i = 0; i < SIZE-1; i++)
	{
		p1 = {pt1[0] + (double)i * dx, pt1[1] + (double)i * dy};
		p2 = {pt1[0] + (double)(i+1) * dx, pt1[1] + (double)(i+1) * dy };
		task.bodies[0].constraint.push_back(make_pair(p1, p2));
	}
	task.bodies[0].constraint.push_back(make_pair(pt1, pt2));
	task.bodies[0].constraint.push_back(make_pair(pt2, pt3));

	dx = (pt4[0] - pt3[0]) / (double)SIZE;
	dy = (pt4[1] - pt3[1]) / (double)SIZE;
	for (int i = 0; i < SIZE-1; i++)
	{
		p1 = { pt3[0] + (double)i * dx, pt3[1] + (double)i * dy };
		p2 = { pt3[0] + (double)(i + 1) * dx, pt3[1] + (double)(i + 1) * dy };
		task.bodies[0].constraint.push_back(make_pair(p1, p2));
	}
	task.bodies[0].constraint.push_back(make_pair(pt3, pt4));
	task.bodies[0].constraint.push_back(make_pair(pt4, pt1));

	oil2d::Oil2d model;
	oil2d::Properties props;
	model.load(task, props);
	model.setSnapshotter(&model);
	//model.snapshotter = std::make_shared<VTKSnapshotter<oil2d::Oil2d>>(&model);
	model.snapshot_all(0);

	return 0;
}
