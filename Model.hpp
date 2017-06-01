#ifndef MODEL_HPP_
#define MODEL_HPP_

#include <memory>

#include "TriangleMesh.hpp"
#include "Variables.hpp"

class FirstModel
{
public:
	typedef mesh::TriangleMesh<var::Var1phase> Mesh;
	std::shared_ptr<Mesh> mesh;
};

#endif /* MODEL_HPP_ */
