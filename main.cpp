#define _USE_MATH_DEFINES
#include <cmath>

#include "src/Scene.hpp"

#include "src/models/Oil2d/Oil2d.hpp"
#include "src/models/Oil2d/Oil2dSolver.hpp"

#include "src/models/Acid/Acid2d.hpp"
#include "src/models/Acid/Acid2dSolver.hpp"

using namespace std;
using namespace point;
using namespace mesh;

/*oil2d::Properties* getProps()
{
	oil2d::Properties* props = new oil2d::Properties();
	props->timePeriods.push_back(86400.0 * 20.0);

	props->leftBoundIsRate = true;
	props->rightBoundIsPres = true;
	props->rates.push_back(1000.0);
	props->ht = 1000.0;
	props->ht_min = 1000.0;
	props->ht_max  = 100000.0;

	props->perfIntervals.push_back( make_pair(0, 0) );
	props->r_w = 0.1;
	props->R_dim = props->r_w * 1.0;
	props->r_e = 1000.0;

	oil2d::Skeleton_Props tmp;
	tmp.m = 0.1;
	tmp.p_init = tmp.p_out = 200.0 * 100000.0;
	tmp.height = 10.0;
	tmp.kx = tmp.ky = 10.0;
	tmp.dens_stc = 2000.0;
	tmp.beta = 0.0 * 1.0e-10;
	props->props_sk.push_back( tmp );

	props->props_oil.beta = 1.e-9;
	props->props_oil.dens_stc = 800.0;
	props->props_oil.visc = 1.0;
	props->props_oil.p_ref = tmp.p_init;

	props->alpha = 7200.0;

	return props;
}*/
acid2d::Properties* getProps()
{
	acid2d::Properties* props = new acid2d::Properties();

	props->timePeriods.push_back(5.0 * 3600.0);
	props->timePeriods.push_back(30.0 * 3600.0);
	props->leftBoundIsRate = true;
	props->rightBoundIsPres = true;
	props->rates.push_back(-10.0);
	props->rates.push_back(0.0);
	//props->pwf.push_back(250.0 * 1.0e+5);
	//props->pwf.push_back(250.0 * 1.0e+5);
	props->xa.push_back(0.13);
	props->xa.push_back(0.0);

	props->ht = 0.1;
	props->ht_min = 0.1;
	props->ht_max = 10000.0;

	props->alpha = 7200.0;

	props->r_w = 0.1;
	props->R_dim = props->r_w;
	props->r_e = 150.0;

	props->perfIntervals.push_back(make_pair(1, 1));
	//props->perfIntervals.push_back(make_pair(15, 15));

	acid2d::Skeleton_Props tmp;
	tmp.cellsNum_z = 1;
	tmp.m_init = 0.1;
	tmp.p_init = tmp.p_out = tmp.p_ref = 200.0 * 1.0e+5;
	tmp.s_init = 0.2;
	tmp.xa_init = 0.0;	tmp.xw_init = 1.0;
	tmp.s_wc = 0.0;		tmp.s_oc = 0.0;		tmp.s_gc = 0.00;
	tmp.xa_eqbm = 0.0;
	tmp.h1 = 0.0;
	tmp.h2 = 10.0;
	tmp.height = 10.0;
	tmp.kx = tmp.ky = 10.0;
	tmp.dens_stc = 2000.0;
	tmp.beta = 4.35113e-10;

	props->props_sk.push_back(tmp);

	props->depth_point = 0.0;

	props->props_o.visc = 1.0;
	props->props_o.dens_stc = 887.261;
	props->props_o.beta = 1.0 * 1.e-9;
	props->props_o.p_ref = tmp.p_ref;

	props->props_w.visc = 1.0;
	props->props_w.dens_stc = 1000.0;
	props->props_w.beta = 1.0 * 1.e-9;
	props->props_w.p_ref = tmp.p_ref;

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

	const double w = 1 / x_dim;
	Point pt1 = { -100 / x_dim, w / 2.0 / x_dim };			Point pt2 = { 100 / x_dim, w / 2.0 / x_dim };
	Point pt3 = pt2;		pt3[1] -= w;
	Point pt4 = pt1;		pt4[1] -= w;

	const int SIZE = 5000;
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

double acid2d::Component::T = 300.0;

int main(int argc, char* argv[])
{
	const auto props = getProps();
	const auto task = getMeshTask(props->R_dim);
	
	//Scene<oil2d::Oil2d, oil2d::Oil2dSolver, oil2d::Properties> scene;
	Scene<acid2d::Acid2d, acid2d::Acid2dSolver, acid2d::Properties> scene;

	scene.load(*props, *task);
	scene.start();

	return 0;
}
