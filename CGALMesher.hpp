#ifndef CGALMESHER_HPP_
#define CGALMESHER_HPP_

#include <string>
#include <fstream>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Polygon_2.h>

namespace cgalmesher
{
	typedef CGAL::Exact_predicates_inexact_constructions_kernel    K;
	typedef CGAL::Triangulation_vertex_base_with_info_2<size_t, K> Vb;
	typedef CGAL::Triangulation_face_base_with_info_2<size_t, K>   Cb;
	typedef CGAL::Triangulation_data_structure_2<Vb, Cb>           Tds;

	/// The result of meshing to be copied to some resulting triangulation
	typedef CGAL::Delaunay_triangulation_2<K, Tds> IntermediateTriangulation;
	typedef IntermediateTriangulation::Point                       CgalPoint2;
	typedef IntermediateTriangulation::Geom_traits::Vector_2       CgalVector2;
	typedef CGAL::Polygon_2<K, std::vector<CgalPoint2>>            Polygon;
	
	/**
	* 2D mesher by CGAL library
	*/
	class Cgal2DMesher {
	public:
		/// Body structure used in mesher task
		struct TaskBody {
			typedef std::array<double, 2> Point;
			/// closed simple polygon without self-intersections
			typedef std::vector<Point>    Border;

			size_t id;                 ///< will be set to cells info
			Border outer;              ///< outer border of the body
			std::vector<Border> inner; ///< borders of the inner cavities of the body
		};

		/// Body structure used inside the mesher
		struct CgalBody {
			size_t id;                  ///< will be set to cells info
			Polygon outer;              ///< outer border of the body
			std::vector<Polygon> inner; ///< borders of the inner cavities of the body

			CgalBody(const TaskBody& task) :
				id(task.id), outer(makePolygon(task.outer)) {
				for (const auto& cavity : task.inner) {
					inner.push_back(makePolygon(cavity));
				}
			}

			bool contains(const CgalPoint2& p) const {
				if (outer.has_on_unbounded_side(p)) { return false; }
				for (const Polygon& cavity : inner) {
					if (cavity.has_on_bounded_side(p)) { return false; }
				}
				return true;
			}
		};


		/// @name Converter structures for CGAL copy_tds function --
		/// convertion between vertices and cells from
		/// different types of triangulations

		template<typename Src, typename Res>
		struct DefaultVertexConverter {
			Res operator()(const Src& src) const {
				return Res(src.point()); // copy coordinates
			}
			void operator()(const Src&, Res&) const { }
		};

		template<typename Src, typename Res>
		struct DefaultCellConverter {
			Res operator()(const Src& src) const {
				Res res;
				// copy cell info
				res.info().setGridId(src.info());
				return res;
			}
			void operator()(const Src&, Res&) const { }
		};

		template<typename Src, typename Res, typename Tr>
		struct BodiesCellConverter {

			Res operator()(const Src& src) const {
				Res res;
				res.info() = (size_t)(-1);

				if (!src.is_in_domain()) { return res; }

				int containsCounter = 0;
				for (const auto& body : *bodies) {
					if (body.contains(centroid(src))) {
						res.info() = body.id;
						++containsCounter;
					}
				}
				assert(containsCounter == 1);

				return res;
			}

			void operator()(const Src&, Res&) const { }

			BodiesCellConverter(const std::vector<CgalBody>* bodies_,
				const Tr* tr_) : bodies(bodies_), tr(tr_) { }

		private:
			const std::vector<CgalBody>* bodies;
			const Tr* const tr;
		};

		/// @}


		/**
		* Build the grid on given geometry
		* @param spatialStep         effective spatial step
		* @param bodies              list of bodies to construct
		* @param result              triangulation to write the result in
		* @tparam ResultingTriangulation type of the triangulation to write result in
		* @tparam CellConverter          see DefaultCellConverter
		* @tparam VertexConverter        see DefaultVertexConverter
		*/
		template<
			typename ResultingTriangulation,
			template<typename, typename> class CellConverter = DefaultCellConverter,
			template<typename, typename> class VertexConverter = DefaultVertexConverter
		>
			static void triangulate(
				const double spatialStep, const std::vector<TaskBody> bodies,
				ResultingTriangulation& result) {

			copyTriangulation<IntermediateTriangulation, ResultingTriangulation,
				CellConverter, VertexConverter>(
					triangulate(spatialStep, convert(bodies)), result);
		}


	private:
		/** The meshing itself */
		static IntermediateTriangulation triangulate(
			const double spatialStep, const std::vector<CgalBody> bodies);


		/**
		* Copy CGAL triangulations of different types
		*/
		template<
			typename InputTriangulation,
			typename OutputTriangulation,
			template<typename, typename> class CellConverter,
			template<typename, typename> class VertexConverter
		>
			static void copyTriangulation(const InputTriangulation& input,
				OutputTriangulation& output) {

			typedef typename InputTriangulation::Face     InputCell;
			typedef typename InputTriangulation::Vertex   InputVertex;
			typedef typename OutputTriangulation::Face    OutputCell;
			typedef typename OutputTriangulation::Vertex  OutputVertex;

			CellConverter<InputCell, OutputCell>        cellConverter;
			VertexConverter<InputVertex, OutputVertex>  vertexConverter;

			copyTriangulation(input, output, cellConverter, vertexConverter);
		}


		/**
		* Copy CGAL triangulations of different types
		*/
		template<
			typename InputTriangulation,
			typename OutputTriangulation,
			typename CellConverter,
			typename VertexConverter
		>
			static void copyTriangulation(
				const InputTriangulation& input, OutputTriangulation& output,
				const CellConverter& cellConverter,
				const VertexConverter& vertexConverter) {
			// try to repeat Triangulation_2 copy constructor as much as possible
			output.set_infinite_vertex(output.tds().copy_tds(
				input.tds(), input.infinite_vertex(), vertexConverter, cellConverter));
		}


		/// @name help functions for mesher
		/// @{
		static Polygon makePolygon(const TaskBody::Border& points);
		static CgalPoint2 findInnerPoint(const Polygon& polygon);

		static std::vector<CgalBody> convert(const std::vector<TaskBody>& bodies) {
			std::vector<CgalBody> ans;
			for (const auto& body : bodies) { ans.push_back(body); }
			return ans;
		}

		/// insert points and constraints - lines between them - to cdt
		/// @tparam CDT Constrained Delaunay triangulation
		template<typename CDT>
		static void insertPolygon(const Polygon& polygon, CDT& cdt) {
			auto point = polygon.vertices_begin();
			typename CDT::Vertex_handle first = cdt.insert(*point);
			typename CDT::Vertex_handle last = first;
			++point;

			while (point != polygon.vertices_end()) {
				typename CDT::Vertex_handle current = cdt.insert(*point);
				cdt.insert_constraint(last, current);
				last = current;
				++point;
			}

			cdt.insert_constraint(last, first);
		}

		/// compute barycenter of the given cell
		template<typename Cell>
		static CgalPoint2 centroid(const Cell& cell) {
			return CGAL::centroid(cell.vertex(0)->point(),
				cell.vertex(1)->point(), cell.vertex(2)->point());
		}

		/// @}
	};
};

#endif /* CGALMESHER_HPP_ */