#ifndef TRIANGLEMESH_HPP_
#define TRIANGLEMESH_HPP_

#include <array>
#include <valarray>
#include <set>
#include <utility>
#include <CGAL/Triangle_2.h>
#include <CGAL/Polygon_2.h>

#include "src/models/Variables.hpp"
#include "src/models/Cell.hpp"
#include "src/models/Element.hpp"
#include "src/mesh/CGALMesher.hpp"
#include "src/snapshotter/VTKSnapshotter.hpp"

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
		template<typename> friend class VTKSnapshotter;
		template<typename> friend class AbstractSolver;
		friend class Oil2dSolver;
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

		typedef TriangleCell Cell;

		static const int CELL_POINTS_NUMBER = 3;	
	protected:
		Triangulation triangulation;
		size_t inner_cells = 0, inner_beg;
		size_t border_edges = 0, border_beg;
		size_t constrained_edges = 0, constrained_beg;
		size_t well_idx;
		std::vector<TriangleCell> cells;
		std::vector<VertexHandle> vertexHandles;
		std::vector<TriangleCell*> fracCells;
		double Volume;

		void load(const Task& task)
		{
			// Task reading
			typedef cgalmesher::Cgal2DMesher::TaskBody Body;
			std::vector<Body> bodies;
			for (const auto& b : task.bodies)
				bodies.push_back({ b.id, b.outer, b.inner, b.constraint });

			// Triangulation sends addtional info about constraints
			std::vector<size_t> constrainedCells;
			cgalmesher::Cgal2DMesher::triangulate(task.spatialStep, bodies, triangulation, constrainedCells);

			// Cells / Vertices addition
			std::set<VertexHandle> localVertices;
			size_t cell_idx = 0;
			inner_beg = 0;
			for (auto cellIter = triangulation.finite_faces_begin(); cellIter != triangulation.finite_faces_end(); ++cellIter)
			{
				cellIter->info().id = cell_idx;
				cells.push_back(TriangleCell(cell_idx++));
				for (int i = 0; i < CELL_POINTS_NUMBER; i++)
					localVertices.insert(cellIter->vertex(i));

				auto& cell = cells[cells.size() - 1];
				cell.type = CellType::INNER;
				const auto& tri = triangulation.triangle(cellIter);
				cell.V = fabs(tri.area());
				Volume += cell.V;
				const auto center = CGAL::barycenter(tri.vertex(0), 1.0 / 3.0, tri.vertex(1), 1.0 / 3.0, tri.vertex(2));
				cell.c = { center[0], center[1] };
			}
			inner_cells = cells.size();

			// Coping vertices from set to vector
			vertexHandles.assign(localVertices.begin(), localVertices.end());
			for (size_t i = 0; i < vertexHandles.size(); i++)
				vertexHandles[i]->info() = i;

			// Border / Constrained cells treatment
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
						const point::Point2d& p1 = { cellIter->vertex(cellIter->cw(i))->point()[0], cellIter->vertex(cellIter->cw(i))->point()[1] };
						const point::Point2d& p2 = { cellIter->vertex(cellIter->ccw(i))->point()[0], cellIter->vertex(cellIter->ccw(i))->point()[1] };
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
				/*if (vecCounter < constrainedCells.size())
				{
					if (constrainedCells[vecCounter] == cell_idx)
					{
						vecCounter++;
						cell.type = CellType::CONSTRAINED;
					}
				}
				if (nebrCounter == CELL_POINTS_NUMBER)
					cell.type = CellType::INNER;
				else if (nebrCounter < CELL_POINTS_NUMBER)
					cell.type = CellType::INNER;*/

				++cellIter;
			}

			// Setting type to fracture cells
			typedef Triangulation::Point CgalPoint2;
			typedef CGAL::Polygon_2<K, std::vector<CgalPoint2>>	Polygon;
			Polygon frac;
			for (const auto& con : task.bodies[0].constraint)
				frac.push_back(CgalPoint2(con.first[0], con.first[1]));
			cell_idx = 0;
			point::Point2d well_pos = { 0.0, 0.0 };
			for (auto cellIter = triangulation.finite_faces_begin(); cellIter != triangulation.finite_faces_end(); ++cellIter)
			{
				auto& cell = cells[cell_idx];
				CgalPoint2 center(cell.c.x, cell.c.y);
				if (!frac.has_on_unbounded_side(center))
				{
					cell.type = CellType::FRAC;
					fracCells.push_back(&cell);
					well_pos.x += center[0];		well_pos.y += center[1];
				}

				cell_idx++;
			}
			well_pos.x /= fracCells.size();		well_pos.y /= fracCells.size();

			std::vector<CellHandle> well_faces;
			std::back_insert_iterator<std::vector<CellHandle>> it(well_faces);
			it = triangulation.get_conflicts(CgalPoint2(well_pos.x, well_pos.y), std::back_inserter(well_faces));

			well_idx = well_faces[0]->info().id;
			cells[well_idx].type = CellType::WELL;
			// Creating constrained cells
/*			constrained_beg = cells.size();
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
			}*/
		};
	public:
		TriangleMesh() { Volume = 0.0; };
		TriangleMesh(const Task& task)
		{
			Volume = 0.0;
			load(task);
		};
		~TriangleMesh() {};

		int getCellsSize() const
		{
			return cells.size();
		}
		int getVerticesSize() const
		{
			return vertexHandles.size();
		}
	};
};

#endif /* TRIANGLEMESH_HPP_ */
