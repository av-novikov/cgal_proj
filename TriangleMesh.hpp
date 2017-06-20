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
		/// The triangulation type
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
		std::vector<VertexHandle> vertexHandles;
		std::vector<CellHandle> cellHandles;
		std::vector<elem::Edge> borderEdges;
		std::vector<CellHandle*> innerHandles;
		std::vector<CellHandle*> borderHandles; 
		std::vector<CellHandle*> constrainedHandles;

		std::shared_ptr<VTKSnapshotter<FirstModel>> snapshotter;

		/*std::list<CellHandle> localIncidentCells(const LocalVertexIndex it) const {
			VertexHandle vh = vertexHandle(it);
			std::list<CellHandle> ans = triangulation->allIncidentCells(vh);
			auto listIter = ans.begin();
			while (listIter != ans.end()) {
				if (belongsToTheGrid(*listIter)) {
					++listIter;
				}
				else {
					listIter = ans.erase(listIter);
				}
			}
			return ans;
		}*/
		/*CellType cellState(const LocalVertexIndex it) const 
		{
			std::set<GridId> incidentGrids = gridsAroundVertex(it);
			assert_true(incidentGrids.erase(id));
			if (incidentGrids.empty()) { return BorderState::INNER; }
			if (incidentGrids.size() > 1) { return BorderState::MULTICONTACT; }
			if (*incidentGrids.begin() == EmptySpaceFlag) { return BorderState::BORDER; }
			return BorderState::CONTACT;
		}*/

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

			std::vector<size_t> constrainedCells;
			cgalmesher::Cgal2DMesher::triangulate(task.spatialStep, bodies, triangulation, constrainedCells);

			std::set<VertexHandle> localVertices;
			size_t i = 0;
			for (auto cellIter = triangulation.finite_faces_begin(); cellIter != triangulation.finite_faces_end(); ++cellIter) 
			{
					cellHandles.push_back(cellIter);
					for (int i = 0; i < CELL_POINTS_NUMBER; i++)
						localVertices.insert(cellIter->vertex(i));
					
					auto& cell = cellHandles[cellHandles.size() - 1];
					cell->info().id = i++;
					const auto& tri = triangulation.triangle(cellIter);
					cell->info().V = fabs(tri.area());
					const auto center = CGAL::barycenter(tri.vertex(0), 1.0 / 3.0, tri.vertex(1), 1.0 / 3.0, tri.vertex(2));
					cell->info().center = { center[0], center[1] };
			}
			vertexHandles.assign(localVertices.begin(), localVertices.end());

			for (size_t i = 0; i < vertexHandles.size(); i++)
				vertexHandles[i]->info() = i;

			int nebrCounter;
			size_t cellCounter = 0;
			size_t vecCounter = 0;
			for (CellIterator cell = cellHandles.begin(); cell != cellHandles.end(); ++cell)
			{
				nebrCounter = 0;
				for (int i = 0; i < CELL_POINTS_NUMBER; i++)
				{
					(*cell)->info().vertexIndices[i] = (*cell)->vertex(i)->info();
					const auto& nebr = (*cell)->neighbor(i);
					if (!triangulation.is_infinite(nebr))
					{
						(*cell)->info().neighborIndices.push_back(nebr->info().id);
						nebrCounter++;
					}
					else
					{
						const Point2d& p1 = { (*cell)->vertex((*cell)->cw(i))->point()[0], (*cell)->vertex((*cell)->cw(i))->point()[1] };
						const Point2d& p2 = { (*cell)->vertex((*cell)->ccw(i))->point()[0], (*cell)->vertex((*cell)->ccw(i))->point()[1] };
						borderEdges.push_back({ point::distance(p1, p2), (p1 + p2) / 2.0,  (*cell)->info().id });
						(*cell)->info().edgesIndices.push_back(borderEdges.size()-1);
					}
				}
				
				if (vecCounter < constrainedCells.size())
				{
					if (constrainedCells[vecCounter] == cellCounter)
					{
						vecCounter++;
						(*cell)->info().type = CellType::CONSTRAINED;
					}
				}

				if(nebrCounter == CELL_POINTS_NUMBER && (*cell)->info().type != CellType::CONSTRAINED)
					(*cell)->info().type = CellType::INNER;
				else if (nebrCounter == CELL_POINTS_NUMBER - 1)
					(*cell)->info().type = CellType::BORDER1;
				else if (nebrCounter == CELL_POINTS_NUMBER - 2)
					(*cell)->info().type = CellType::BORDER2;
				else if (nebrCounter == CELL_POINTS_NUMBER - 3)
					(*cell)->info().type = CellType::BORDER3;

				cellCounter++;
			}
		};
		void snapshot(const int i) const
		{
			snapshotter->dump(i);
		}
	};
};

#endif /* TRIANGLEMESH_HPP_ */
