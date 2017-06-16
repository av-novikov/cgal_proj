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
	task.spatialStep = 0.2;
	Task::Body::Border body1border = { { 3, 3 },{ -3, 3 },{ -3, -3 },{ 3, -3 } };
	Task::Body::Border body2border = { { 3, 3 },{ 9, 3 },{ 9, -3 },{ 3, -3 } };
	task.bodies = {	Task::Body({ 0, body1border, {}})};

	Point pt1 = { 1.0, -1.0 };			Point pt2 = { -1.0, 1.0 };

	Task::Body::Edge frac = make_pair(pt1, pt2);
	task.bodies[0].constraint.push_back(frac);
	const double alpha1 = atan((pt2[1] - pt1[1]) / (pt2[0] - pt1[0])) - M_PI / 2.0;
	const double alpha2 = atan((pt2[1] - pt1[1]) / (pt2[0] - pt1[0])) + M_PI / 2.0;
	const int SIZE = 50;
	const double dalpha = (alpha2 - alpha1) / (double)SIZE;
	const double length = sqrt((pt2[0] - pt1[0]) * (pt2[0] - pt1[0]) + (pt2[1] - pt1[1]) * (pt2[1] - pt1[1]));
	double rad = length / 4.0;
	Point p1;
	task.bodies[0].inner.resize(5);
	for (int k = 0; k < 5; k++)
	{
		rad = length / pow(2.0, (double)k + 2.0);
		for (int i = 0; i < SIZE + 1; i++)
		{
			p1 = { pt1[0] + rad * cos(alpha1 + (double)i * dalpha), pt1[1] + rad * sin(alpha1 + (double)i * dalpha) };
			task.bodies[0].inner[k].push_back(p1);
		}
		for (int i = 0; i < SIZE + 1; i++)
		{
			p1 = { pt2[0] + rad * cos(alpha2 + (double)i * dalpha), pt2[1] + rad * sin(alpha2 + (double)i * dalpha) };
			task.bodies[0].inner[k].push_back(p1);
		}
	}

	FirstModel model;
	model.mesh = make_shared<FirstModel::Mesh>(task);
	model.mesh->snapshot(0);

	return 0;
}
