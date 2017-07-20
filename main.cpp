#define _USE_MATH_DEFINES
#include <cmath>

#include "src/util/utils.h"
#include "src/Scene.hpp"

#include "src/models/Oil2d/Oil2d.hpp"
#include "src/models/Oil2d/Oil2dSolver.hpp"

#include "src/models/Acid/Acid2d.hpp"
#include "src/models/Acid/Acid2dSolver.hpp"

using namespace std;
using namespace point;
using namespace mesh;

typedef Task::Body::Point Point;

#define CIRC_NUM 20

/*oil2d::Properties* getProps()
{
	oil2d::Properties* props = new oil2d::Properties();
	props->timePeriods.push_back(86400.0 * 20.0);

	props->leftBoundIsRate = true;
	props->rightBoundIsPres = true;
	props->rates.push_back(100.0);
	props->ht = 100.0;
	props->ht_min = 100.0;
	props->ht_max  = 100000.0;

	props->perfIntervals.push_back( make_pair(0, 0) );
	props->r_w = 0.1;
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
	props->rates.push_back(10.0);
	props->rates.push_back(0.0);
	//props->pwf.push_back(250.0 * 1.0e+5);
	//props->pwf.push_back(250.0 * 1.0e+5);
	props->xa.push_back(0.13);
	props->xa.push_back(0.0);

	props->ht = 10.0;
	props->ht_min = 10.0;
	props->ht_max = 10000.0;

	props->alpha = 7200.0;

	props->r_w = 0.1;
	props->R_dim = props->r_w * 10.0;
	props->r_e = 150.0;

	props->perfIntervals.push_back(make_pair(1, 1));
	//props->perfIntervals.push_back(make_pair(15, 15));

	acid2d::Skeleton_Props tmp;
	tmp.cellsNum_z = 1;
	tmp.m_init = 0.1;
	tmp.p_init = tmp.p_out = tmp.p_ref = 200.0 * 1.0e+5;
	tmp.s_init = 0.8;
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

void getFracPoints(Point2d point1, Point2d point2, const double w, const double r_w, vector<Point>& pts)
{
	const size_t size = pts.size();
	pts.clear();
	const double k = (point2.y - point1.y) / (point2.x - point1.x);
	const double alpha = atan(k);
	const double l = distance(point1, point2) / 2.0;

	double z_prev1 = l * pow(l / r_w, -2.0 / (double)(size - 1)), z_prev2 = r_w;
	double logMax_z = log(l / r_w);
	double logStep_z = 2.0 * logMax_z / (double)(size - 1);

	Point pt1 = { point1.x + w / 2.0 * sin(alpha), point1.y - w / 2.0 * cos(alpha) };
	Point pt2 = { point2.x + w / 2.0 * sin(alpha), point2.y - w / 2.0 * cos(alpha) };
	Point pt3 = pt2;		pt3[0] -= w * sin(alpha);		pt3[1] += w * cos(alpha);
	Point pt4 = pt1;		pt4[0] -= w * sin(alpha);		pt4[1] += w * cos(alpha);
	Point2d well_pt = (point1 + point2) / 2.0;

	Point p1;
	double upcoord, hz = 0.0, cm_z = 0.0, hz_prev = 0.0;
	p1 = pt1;
	pts.push_back({ pt1[0], pt1[1] });
	for (int i = 0; i < size; i++)
	{
		upcoord = cm_z + hz / 2.0 + r_w;
		if (upcoord < l - EQUALITY_TOLERANCE)
		{
			hz = z_prev1 * (exp(logStep_z) - 1.0);
			z_prev1 *= exp(-logStep_z);
		}
		else if (upcoord < l + r_w - EQUALITY_TOLERANCE)
		{
			hz = 2.0 * r_w;
			const double hphi = (pts[pts.size()-1][0] - pts[pts.size() - 2][0] > 0.0) ? M_PI / (double)(CIRC_NUM + 1) : -M_PI / (double)(CIRC_NUM + 1);
			const double start_phi = (point2.x - point1.x < 0.0) ? alpha : M_PI + alpha;
			for (int j = 1; j <= CIRC_NUM; j++)
				pts.push_back({ well_pt.x + r_w * cos(start_phi + (double)j * hphi), well_pt.y + r_w * sin(start_phi + (double)j * hphi) });
		}
		else
		{
			hz = z_prev2 * (exp(logStep_z) - 1.0);
			z_prev2 *= exp(logStep_z);
		}

		cm_z += (hz_prev + hz) / 2.0;
		hz_prev = hz;

		pts.push_back({ pt1[0] + sign(point2.x - point1.x) * (cm_z + hz / 2.0) * cos(alpha), pt1[1] + sign(point2.y - point1.y) * (cm_z + hz / 2.0) * sin(fabs(alpha)) });
		//p2 = { pt1[0] + (cm_z + hz / 2.0) * cos(alpha), pt1[1] + (cm_z + hz / 2.0) * sin(alpha) };
		//task->bodies[0].constraint.push_back(make_pair(p1, p2));
		//p1 = p2;
	}
	pts.push_back({ pt3[0], pt3[1] });
	//task->bodies[0].constraint.push_back(make_pair(p2, pt3));

	cm_z = 0.0;
	p1 = pt3;
	z_prev1 = l * pow(l / r_w, -2.0 / (double)(size - 1));
	z_prev2 = r_w;
	hz_prev = hz = 0.0;
	for (int i = 0; i < size; i++)
	{
		upcoord = cm_z + hz / 2.0 + r_w;
		if (upcoord < l - EQUALITY_TOLERANCE)
		{
			hz = z_prev1 * (exp(logStep_z) - 1.0);
			z_prev1 *= exp(-logStep_z);
		}
		else if (upcoord < l + r_w - EQUALITY_TOLERANCE)
		{
			hz = 2.0 * r_w;
			const double hphi = (pts[pts.size() - 2][0] - pts[pts.size() - 1][0] > 0.0) ? M_PI / (double)(CIRC_NUM + 1) : -M_PI / (double)(CIRC_NUM + 1);
			const double start_phi = (point1.x - point2.x < 0.0) ? alpha : M_PI + alpha;
			for (int j = 1; j <= CIRC_NUM; j++)
				pts.push_back({ well_pt.x + r_w * cos(start_phi + (double)j * hphi), well_pt.y + r_w * sin(start_phi + (double)j * hphi) });
		}
		else
		{
			hz = z_prev2 * (exp(logStep_z) - 1.0);
			z_prev2 *= exp(logStep_z);
		}

		cm_z += (hz_prev + hz) / 2.0;
		hz_prev = hz;

		pts.push_back({ pt3[0] - sign(point2.x - point1.x) * (cm_z + hz / 2.0) * cos(alpha), pt3[1] - sign(point2.y - point1.y) * (cm_z + hz / 2.0) * sin(fabs(alpha)) });
		//p2 = { pt3[0] - (cm_z + hz / 2.0) * cos(alpha), pt3[1] - (cm_z + hz / 2.0) * sin(alpha) };
		//task->bodies[0].constraint.push_back(make_pair(p1, p2));
		//p1 = p2;
	}
	//task->bodies[0].constraint.push_back(make_pair(p2, pt1));
}
Task* getMeshTask(double& x_dim, double r_w)
{
	Task* task = new Task;

	double w = 0.01;
	x_dim = w * 20;		w /= x_dim;
	r_w /= x_dim;
	const double a = 300.0 / x_dim;
	const double l = 100.0 / x_dim;

	task->spatialStep = 50.0 / x_dim;
	Point2d pt1 = { l, l }, pt2 = { -l, -l };
	Point2d pt_well = (pt1 + pt2) / 2.0;
	Task::Body::Border body1border = { { a, a }, { -a, a }, { -a, -a }, { a, -a } };
	task->bodies = { Task::Body({ 0, r_w, {pt_well.x, pt_well.y}, body1border,{} }) };

	const size_t size = 101;
	vector<Point> pts(size);
	getFracPoints(pt1, pt2, w, r_w, pts);
	for(size_t i = 1; i < pts.size(); i++)
		task->bodies[0].constraint.push_back(make_pair(pts[i-1], pts[i]));
	task->bodies[0].constraint.push_back(make_pair(pts[pts.size()-1], pts[0]));

	return task;
}

double acid2d::Component::T = 300.0;

int main(int argc, char* argv[])
{
	const auto props = getProps();
	const auto task = getMeshTask(props->R_dim, props->r_w);
	
	//Scene<oil2d::Oil2d, oil2d::Oil2dSolver, oil2d::Properties> scene;
	Scene<acid2d::Acid2d, acid2d::Acid2dSolver, acid2d::Properties> scene;

	scene.load(*props, *task);
	scene.start();

	return 0;
}
