#ifndef TRIANGLEMESH_HPP_
#define TRIANGLEMESH_HPP_

#include "Cell.hpp"
#include "AbstractMesh.hpp"

namespace mesh
{
	template <typename TVariable>
	class TriangleMesh : public AbstractMesh<cell::TriangleCell<TVariable> >
	{

	};
};

#endif /* TRIANGLEMESH_HPP_ */
