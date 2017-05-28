#ifndef ABSTRACTMESH_HPP_
#define ABSTRACTMESH_HPP_

#include <vector>
#include <valarray>

namespace mesh
{
	template <class TCell>
	class AbstractMesh
	{
	public:
		typedef TCell Cell;
		typedef typename Cell::Point Point;
		typedef typename Cell::Element Element;
		typedef typename Cell::Facet Facet;
	protected:
		std::valarray<Cell> cells;
		std::valarray<Point> points;
		std::valarray<Facet> facets;
		// Slices !!!
		std::vector<Facet*> borderFacets;
		std::vector<Facet*> sourceFacets;
	public:

		AbstractMesh() {};
		virtual ~AbstractMesh() {};
	};
};

#endif /* ABSTRACTMESH_HPP_ */