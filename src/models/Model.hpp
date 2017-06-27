#ifndef MODEL_HPP_
#define MODEL_HPP_

#include <memory>

#include "src/mesh/TriangleMesh.hpp"
#include "src/models/Variables.hpp"

class FirstModel
{
public:
	typedef var::containers::Var1phase Variable;
	typedef mesh::TriangleMesh<Variable> Mesh;
	std::shared_ptr<Mesh> mesh;
};

#endif /* MODEL_HPP_ */
