#ifndef TRIANGLEMESH_HPP_
#define TRIANGLEMESH_HPP_

#include <array>
#include <set>

#include "Cell.hpp"
#include "AbstractMesh.hpp"
#include "CGALMesher.hpp"
#include "VTKSnapshotter.hpp"

struct Task
{
	double spatialStep;
	struct Body {
		typedef std::array<double, 2> Point;
		typedef std::vector<Point> Border;

		size_t id;                  ///< body indicator > 0 @see Task::Body
		Border outer;               ///< outer border of the body
		std::vector<Border> inner;  ///< borders of the inner cavities
	};
	std::vector<Body> bodies;
};

class FirstModel;

namespace mesh
{
	enum CellType { INNER, BORDER_IN, BORDER_OUT };
	struct CellInfo
	{
		size_t localVertexIndices[3];
		size_t localNeighborIndices[3];
		CellType type;
	};
	typedef size_t VertexInfo;
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
	class TriangleMesh : public AbstractMesh<cell::TriangleCell<TVariable> >
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
		std::vector<CellHandle> cellHandles;
		std::vector<VertexHandle> vertexHandles;
		std::vector<LocalVertexIndex> borderIndices;
		std::vector<LocalVertexIndex> innerIndices;

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
				bodies.push_back({ b.id, b.outer, b.inner });

			cgalmesher::Cgal2DMesher::triangulate(task.spatialStep, bodies, triangulation);

			std::set<VertexHandle> localVertices;
			for (auto cellIter = triangulation.finite_faces_begin(); cellIter != triangulation.finite_faces_end(); ++cellIter) 
			{
				//if (cellIter->info().getGridId() == id) 
				//{
					cellHandles.push_back(cellIter);
					for (int i = 0; i < CELL_POINTS_NUMBER; i++)
						localVertices.insert(cellIter->vertex(i));
				//}
			}
			vertexHandles.assign(localVertices.begin(), localVertices.end());

			for (size_t i = 0; i < vertexHandles.size(); i++)
				vertexHandles[i]->info() = i;
			for (CellIterator cell = cellHandles.begin(); cell != cellHandles.end(); ++cell)
				for (int i = 0; i < CELL_POINTS_NUMBER; i++)
					(*cell)->info().localVertexIndices[i] = (*cell)->vertex(i)->info();


		};
		void snapshot(const int i) const
		{
			snapshotter->dump(i);
		}
	};
};

#endif /* TRIANGLEMESH_HPP_ */
