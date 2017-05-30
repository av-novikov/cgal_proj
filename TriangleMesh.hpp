#ifndef TRIANGLEMESH_HPP_
#define TRIANGLEMESH_HPP_

#include "Cell.hpp"
#include "AbstractMesh.hpp"
#include "CGALMesher.hpp"

struct Task
{
	double spatialStep;
	struct Body {
		typedef point::Point2d Point;
		typedef std::vector<Point> Border;

		size_t id;                  ///< body indicator > 0 @see Task::Body
		Border outer;               ///< outer border of the body
		std::vector<Border> inner;  ///< borders of the inner cavities
	};
	std::vector<Body> bodies;
};

template<int CELL_POINTS_NUMBER>
struct CellInfoT {
	/// Indicator that no grid owns the cell (auxiliary empty cell)
	static const GridId EmptySpaceFlag = (GridId)(-1);

	/// global triangulation cell can belongs to the only one grid
	GridId gridId;
	void setGridId(const GridId gridId_) { gridId = gridId_; }
	GridId getGridId() const { return gridId; }

	/// local indices of the cell's vertices in the order
	/// the same with their pointers (VertexHandles)
	typedef size_t LocalVertexIndex;
	LocalVertexIndex localVertexIndices[CELL_POINTS_NUMBER];
};

/** Auxiliary information stored in global triangulation vertices */
typedef size_t VertexInfo;

namespace mesh
{
	template <typename TVariable>
	class TriangleMesh : public AbstractMesh<cell::TriangleCell<TVariable> >
	{
	public:
		typedef CGAL::Exact_predicates_inexact_constructions_kernel        K;
		typedef CGAL::Triangulation_vertex_base_with_info_2<VertexInfo, K> Vb;
		typedef CGAL::Triangulation_face_base_with_info_2<CellInfo, K>     Cb;
		typedef CGAL::Triangulation_data_structure_2<Vb, Cb>               Tds;
		/// The triangulation type
		typedef CGAL::Delaunay_triangulation_2<K, Tds>                     Triangulation;

		typedef typename Triangulation::Vertex_handle           VertexHandle;
		typedef typename Triangulation::Face_handle             CellHandle;
		typedef typename Triangulation::Geom_traits::Vector_2   CgalVectorD;
		typedef typename Triangulation::Point                   CgalPointD;
		typedef typename Triangulation::All_faces_iterator      AllCellsIterator;
		typedef typename Triangulation::Line_face_circulator    LineFaceCirculator;

	public:
		TriangleMesh() {};
		TriangleMesh(const Task& task) 
		{
			typedef cgalmesher::Cgal2DMesher::TaskBody Body;
			std::vector<Body> bodies;
			for (const auto& b : task.bodies)
				bodies.push_back({ b.id, b.outer, b.inner });

			cgalmesher::Cgal2DMesher::triangulate(task.spatialStep, bodies, triangulation);
		};
		~TriangleMesh() {};

		Triangulation triangulation;
	};
};

#endif /* TRIANGLEMESH_HPP_ */
