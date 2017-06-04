#include "Model.hpp"

using namespace std;
using namespace point;
using namespace mesh;

int main(int argc, char* argv[])
{
	Task task;
	task.spatialStep = 0.5;
	Task::Body::Border body1border = { { 3, 3 },{ -3, 3 },{ -3, -3 },{ 3, -3 } };
	Task::Body::Border body2border = { { 3, 3 },{ 9, 3 },{ 9, -3 },{ 3, -3 } };
	task.bodies = {	Task::Body({ 0, body1border,{} })};

	FirstModel model;
	model.mesh = make_shared<FirstModel::Mesh>(task);
	model.mesh->snapshot(0);

	return 0;
}
