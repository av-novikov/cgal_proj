#ifndef TRIANGLEMESH_HPP_
#define TRIANGLEMESH_HPP_

#include <array>
#include <set>
#include <CGAL/Triangle_2.h>

#include "Cell.hpp"
#include "Element.hpp"
#include "AbstractMesh.hpp"
#include "CGALMesher.hpp"
#include "VTKSnapshotter.hpp"
#include <utility>

class FirstModel;

struct Task
{
	double spatialStep;
	struct Body {
		typedef std::array<double, 2> Point;
		typedef std::vector<Point> Border;
		typedef std::pair<Point, Point> Edge;

		size_t id;                  ///< body indicator > 0 @see Task::Body
		Border outer;               ///< outer border of the body
		std::vector<Border> inner;  ///< borders of the inner cavities
		std::vector<Edge> constraint;
	};
	std::vector<Body> bodies;
};

struct CellInfo
{
	size_t id;
	point::Point2d center;
	double V;
	size_t vertexIndices[3];
	std::vector<size_t> neighborIndices;
	std::vector<size_t> edgesIndices;		// only for border cells
	CellType type;
	bool isConstrained = false;
};
typedef size_t VertexInfo;

namespace mesh
{
	struct Iterator {
		typedef size_t Index;

		Index iter = 0;

		Iterator(size_t value = 0) : iter(value) { }

		operator Index() const { return iter; }

		const Iterator& operator*() const { return *this; }

		bool operator==(const Iterator& other) const {
			return iter == other.iter;
		}

		bool operator!=(const Iterator& other) const {
			return !((*this) == other);
		}

		bool operator<(const Iterator& other) const {
			return iter < other.iter;
		}

		Iterator& operator++() {
			iter++;
			return (*this);
		}
	};

	template <typename TVariable>
	class TriangleMesh
	{
	public: 
		typedef CGAL::Exact_predicates_inexact_constructions_kernel        K;
		typedef CGAL::Triangulation_vertex_base_with_info_2<VertexInfo, K> Vb;
		typedef CGAL::Triangulation_face_base_with_info_2<CellInfo, K>     Cb;
		typedef CGAL::Triangulation_data_structure_2<Vb, Cb>               Tds;
		typedef CGAL::Delaunay_triangulation_2<K, Tds>                     Triangulation;

		typedef Triangulation::Vertex_handle           VertexHandle;
		typedef Triangulation::Face_handle             CellHandle;
		typedef Triangulation::Geom_traits::Vector_2   CgalVectorD;
		typedef Triangulation::Point                   CgalPointD;
		typedef Triangulation::All_faces_iterator      AllCellsIterator;
		typedef Triangulation::Line_face_circulator    LineFaceCirculator;
		typedef Iterator::Index LocalVertexIndex;
		typedef std::vector<CellHandle>::const_iterator CellIterator;

		static const int CELL_POINTS_NUMBER = 3;	
	public:
		Triangulation triangulation;
		int inner_cells = 0, inner_beg;
		int border_edges = 0, border_beg;
		int constrained_edges = 0, constrained_beg;
		std::vector<TriangleCell> cells;
		std::vector<VertexHandle> vertexHandles;
		std::shared_ptr<VTKSnapshotter<FirstModel>> snapshotter;
	public:
		TriangleMesh() {};
		TriangleMesh(const Task& task) 
		{
			load(task);
			snapshotter = std::make_shared<VTKSnapshotter<FirstModel>>(this);
		};
		~TriangleMesh() {};

		void load(const Task& task)
		{
			typedef cgalmesher::Cgal2DMesher::TaskBody Body;
			std::vector<Body> bodies;
			for (const auto& b : task.bodies)
				bodies.push_back({ b.id, b.outer, b.inner, b.constraint });

			std::vector<std::pair<size_t, int>> constrainedCells;
			cgalmesher::Cgal2DMesher::triangulate(task.spatialStep, bodies, triangulation, constrainedCells);

			std::set<VertexHandle> localVertices;
			size_t cell_idx = 0;
			inner_beg = 0;
			for (auto cellIter = triangulation.finite_faces_begin(); cellIter != triangulation.finite_faces_end(); ++cellIter) 
			{
				cellIter->info().id = cell_idx;
				cells.push_back(TriangleCell(cell_idx++));
				for (int i = 0; i < CELL_POINTS_NUMBER; i++)
					localVertices.insert( cellIter->vertex(i) );
				
				auto& cell = cells[cells.size() - 1];
				cell.type = CellType::INNER;
				const auto& tri = triangulation.triangle(cellIter);
				cell.V = fabs(tri.area());
				const auto center = CGAL::barycenter(tri.vertex(0), 1.0 / 3.0, tri.vertex(1), 1.0 / 3.0, tri.vertex(2));
				cell.c = { center[0], center[1] };
			}
			inner_cells = cells.size();

			vertexHandles.assign(localVertices.begin(), localVertices.end());
			for (size_t i = 0; i < vertexHandles.size(); i++)
				vertexHandles[i]->info() = i;

			//setConstrainedEdges(constrainedEdges);

			int nebrCounter;
			size_t vecCounter = 0;
			auto cellIter = triangulation.finite_faces_begin();
			cells.reserve(int(1.5  * inner_cells));
			border_beg = cells.size();
			for (int cell_idx = 0; cell_idx < inner_cells; cell_idx++)
			{
				TriangleCell& cell = cells[cell_idx];
				nebrCounter = 0;
				for (int i = 0; i < CELL_POINTS_NUMBER; i++)
				{
					cell.points[i] = cellIter->vertex(i)->info();
					const auto& nebr = cellIter->neighbor(i);
					if (!triangulation.is_infinite(nebr))
					{
						cell.nebr[i] = nebr->info().id;
						nebrCounter++;
					}
					else
					{
						const Point2d& p1 = { cellIter->vertex(cellIter->cw(i))->point()[0], cellIter->vertex(cellIter->cw(i))->point()[1] };
						const Point2d& p2 = { cellIter->vertex(cellIter->ccw(i))->point()[0], cellIter->vertex(cellIter->ccw(i))->point()[1] };
						cells.push_back(TriangleCell(inner_cells + border_edges));
						TriangleCell& edge = cells[cells.size() - 1];
						edge.type = CellType::BORDER;
						edge.c = (p1 + p2) / 2.0;
						edge.V = point::distance(p1, p2);
						edge.nebr[0] = cell_idx;
						edge.points[0] = cellIter->vertex(cellIter->cw(i))->info();
						edge.points[1] = cellIter->vertex(cellIter->ccw(i))->info();
						cell.nebr[i] = inner_cells + border_edges;
						border_edges++;
					}
				}
				if (vecCounter < constrainedCells.size())
				{
					if (constrainedCells[vecCounter].first == cell_idx)
					{
						vecCounter++;
						cell.type = CellType::CONSTRAINED;
					}
				}

				if (nebrCounter == CELL_POINTS_NUMBER && cell.type != CellType::CONSTRAINED)
					cell.type = CellType::INNER;
				else if (nebrCounter < CELL_POINTS_NUMBER)
					cell.type = CellType::INNER;

				++cellIter;
			}

			constrained_beg = cells.size();
			cellIter = triangulation.finite_faces_begin();
			for (int cell_idx = 0; cell_idx < inner_cells; cell_idx++)
			{
				auto& cell = cells[cell_idx];
				if (cell.type == CellType::CONSTRAINED)
				{
					const int nebr_idx = constrainedCells[constrained_edges].second;
					const auto nebrIter = cellIter->neighbor(nebr_idx);
					const int nebr_id = nebrIter->info().id;
					auto& nebr = cells[nebrIter->info().id];
					cell.type = nebr.type = CellType::INNER;

					cells.push_back(TriangleCell(constrained_beg + constrained_edges));
					TriangleCell& edge = cells[cells.size() - 1];
					edge.type = CellType::CONSTRAINED;
					const Point2d& p1 = { cellIter->vertex(cellIter->cw(nebr_idx))->point()[0], cellIter->vertex(cellIter->cw(nebr_idx))->point()[1] };
					const Point2d& p2 = { cellIter->vertex(cellIter->ccw(nebr_idx))->point()[0], cellIter->vertex(cellIter->ccw(nebr_idx))->point()[1] };
					edge.c = (p1 + p2) / 2.0;
					edge.V = point::distance(p1, p2);
					edge.nebr[0] = cell_idx;		edge.nebr[1] = nebr.id;
					edge.points[0] = cellIter->vertex(cellIter->cw(nebr_idx))->info();
					edge.points[1] = cellIter->vertex(cellIter->ccw(nebr_idx))->info();
					
					cell.nebr[nebr_idx] = constrained_beg + constrained_edges;
					for (int i = 0; i < CELL_POINTS_NUMBER; i++)
						if (nebr.nebr[i] == cell.id)
							nebr.nebr[i] = cell.nebr[nebr_idx];

					constrained_edges++;
				}

				++cellIter;
			}
		};
		void snapshot(const int i) const
		{
			snapshotter->dump(i);
		}
	};
};

#endif /* TRIANGLEMESH_HPP_ */
