#include "TriangleMesh.hpp"
#include "Variables.hpp"

using namespace std;
using namespace point;
using namespace mesh;

int main(int argc, char* argv[])
{
	Task task;
	task.spatialStep = 0.2;
	Task::Body::Border body1border = { { 3, 3 },{ -3, 3 },{ -3, -3 },{ 3, -3 } };
	Task::Body::Border body2border = { { 3, 3 },{ 9, 3 },{ 9, -3 },{ 3, -3 } };
	task.bodies = {	Task::Body({ 0, body1border,{} }), Task::Body({ 1, body2border,{} }) };

	TriangleMesh<var::Var1phase> mesh(task);

	return 0;
}
