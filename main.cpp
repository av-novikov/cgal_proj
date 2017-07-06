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

	props->leftBoundIsRate = true;
	props->rightBoundIsPres = true;
	props->rates.push_back(100.0);
	props->ht = 10.0;
	props->ht_min = 10.0;
	props->ht_max  = 100000.0;

	props->perfIntervals.push_back( make_pair(0, 0) );
	props->r_w = 0.1;
	props->R_dim = props->r_w * 1.0;
	props->r_e = 1000.0;

	oil2d::Skeleton_Props tmp;
	tmp.m = 0.1;
	tmp.p_init = tmp.p_out = 200.0 * 100000.0;
	tmp.height = 10.0;
	tmp.kx = tmp.ky = 50.0;
	tmp.dens_stc = 2000.0;
	tmp.beta = 0.0 * 1.0e-10;
	props->props_sk.push_back( tmp );

	props->props_oil.beta = 1.e-9;
	props->props_oil.dens_stc = 800.0;
	props->props_oil.visc = 1.0;
	props->props_oil.p_ref = tmp.p_init;

	props->alpha = 7200.0;

	return props;
}

Task* getMeshTask(const double x_dim)
{
	Task* task = new Task;

	typedef Task::Body::Point Point;

	task->spatialStep = 50.0 / x_dim;
	Task::Body::Border body1border = { { 300 / x_dim, 300 / x_dim },
										{ -300 / x_dim, 300 / x_dim },
										{ -300 / x_dim, -300 / x_dim },
										{ 300 / x_dim, -300 / x_dim } };
	//Task::Body::Border body2border = { { 3, 3 },{ 9, 3 },{ 9, -3 },{ 3, -3 } };
	task->bodies = { Task::Body({ 0, body1border,{} }) };

	const double w = 0.1 / x_dim;
	Point pt1 = { -50 / x_dim, -100 / x_dim };			Point pt2 = { -150 / x_dim, 50 / x_dim };
	Point pt3 = pt2;		pt3[0] += w;
	Point pt4 = pt1;		pt4[0] += w;

	const int SIZE = 1;
	double dx = (pt2[0] - pt1[0]) / (double)SIZE;
	double dy = (pt2[1] - pt1[1]) / (double)SIZE;
	Point p1, p2;
	for (int i = 0; i < SIZE - 1; i++)
	{
		p1 = { pt1[0] + (double)i * dx, pt1[1] + (double)i * dy };
		p2 = { pt1[0] + (double)(i + 1) * dx, pt1[1] + (double)(i + 1) * dy };
		task->bodies[0].constraint.push_back(make_pair(p1, p2));
	}
	task->bodies[0].constraint.push_back(make_pair(pt1, pt2));
	task->bodies[0].constraint.push_back(make_pair(pt2, pt3));

	dx = (pt4[0] - pt3[0]) / (double)SIZE;
	dy = (pt4[1] - pt3[1]) / (double)SIZE;
	for (int i = 0; i < SIZE - 1; i++)
	{
		p1 = { pt3[0] + (double)i * dx, pt3[1] + (double)i * dy };
		p2 = { pt3[0] + (double)(i + 1) * dx, pt3[1] + (double)(i + 1) * dy };
		task->bodies[0].constraint.push_back(make_pair(p1, p2));
	}
	task->bodies[0].constraint.push_back(make_pair(pt3, pt4));
	task->bodies[0].constraint.push_back(make_pair(pt4, pt1));

	return task;
}

int main(int argc, char* argv[])
{
	const auto props = getProps();
	const auto task = getMeshTask(props->R_dim);
	Scene<oil2d::Oil2d, oil2d::Oil2dSolver, oil2d::Properties> scene;
	scene.load(*props, *task);
	scene.start();

	return 0;
}
