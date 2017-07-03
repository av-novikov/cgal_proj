#define _USE_MATH_DEFINES
#include <cmath>

#include "src/Scene.hpp"

#include "src/models/Oil2d/Oil2d.hpp"
#include "src/models/Oil2d/Oil2dSolver.hpp"

using namespace std;
using namespace point;
using namespace mesh;

oil2d::Properties* getProps()
{
	oil2d::Properties* props = new oil2d::Properties();
	props->timePeriods.push_back(86400.0 * 20.0);

	props->leftBoundIsRate = false;
	props->rightBoundIsPres = true;
	props->pwf.push_back(180.0 * 1.E+5);
	props->ht = 100.0;
	props->ht_min = 100.0;
	props->ht_max  = 100000.0;

	props->perfIntervals.push_back( make_pair(0, 0) );
	props->r_w = 0.05;
	props->r_e = 1000.0;

	oil2d::Skeleton_Props tmp;
	tmp.m = 0.1;
	tmp.p_init = tmp.p_out = 200.0 * 100000.0;
	tmp.height = 0.1;
	tmp.kx = tmp.ky = 50.0;
	tmp.dens_stc = 2000.0;
	tmp.beta = 0.0 * 1.0e-10;
	props->props_sk.push_back( tmp );

	props->props_oil.beta = 1.e-9;
	props->props_oil.dens_stc = 800.0;
	props->props_oil.visc = 1.0;

	props->alpha = 7200.0;

	return props;
}

int main(int argc, char* argv[])
{
	typedef Task::Body::Point Point;

	Task task;
	task.spatialStep = 50.0;
	Task::Body::Border body1border = { { 300, 300 },{ -300, 300 },{ -300, -300 },{ 300, -300 } };
	//Task::Body::Border body2border = { { 3, 3 },{ 9, 3 },{ 9, -3 },{ 3, -3 } };
	task.bodies = {	Task::Body({ 0, body1border, {}})};

	const double w = 1;
	Point pt1 = { -50, -100 };			Point pt2 = { -150, 50 };
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

	const auto props = getProps();
	Scene<oil2d::Oil2d, oil2d::Oil2dSolver, oil2d::Properties> scene;
	scene.load(*props, task);

	return 0;
}
