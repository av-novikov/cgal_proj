#include <CGAL/Delaunay_mesh_face_base_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Delaunay_mesh_size_criteria_2.h>
#include <CGAL/Delaunay_mesher_2.h>

#include <CGAL/lloyd_optimize_mesh_2.h>
#include "src/mesh/CGALMesher.hpp"

using namespace cgalmesher;

IntermediateTriangulation Cgal2DMesher::triangulate(const double spatialStep, const std::vector<CgalBody> bodies,		
									std::vector<size_t>& constrainedCells)
{
	typedef CGAL::Delaunay_mesh_face_base_2<K>				     Fb;
	typedef CGAL::Triangulation_data_structure_2<Vb, Fb>        Tds;
	typedef CGAL::Constrained_Delaunay_triangulation_2<K, Tds>  CDT;
	typedef CGAL::Delaunay_mesh_size_criteria_2<CDT>            Criteria;
	typedef CGAL::Delaunay_mesher_2<CDT, Criteria>              Mesher;

	CDT cdt;
	// Special "seeds" in the interior of inner cavities
	// to tell CGAL do not mesh these cavities
	std::list<CgalPoint2> listOfSeeds;

	// insert all bodies to triangulation as its constraints
	for (const auto& body : bodies) 
	{
		insertPolygon(body.outer, cdt);
		for (const auto& innerCavity : body.inner) 
		{
			insertPolygon(innerCavity, cdt);
			//listOfSeeds.push_back(findInnerPoint(innerCavity));
		}
		//listOfSeeds.push_back(CDT::Point(1.0, 1.0));
		for (const auto& con : body.constraint)
			cdt.insert_constraint(con.first, con.second);
	}

	Mesher mesher(cdt);
	mesher.set_criteria(Criteria(0.1, spatialStep));
	mesher.refine_mesh();
	//CGAL::lloyd_optimize_mesh_2(cdt, CGAL::parameters::max_iteration_number = 10);

	size_t cell_idx = 0;
	int vert_bounded;
	for (auto cellIter = cdt.finite_faces_begin(); cellIter != cdt.finite_faces_end(); ++cellIter)
	{
		vert_bounded = 0;
		for (int i = 0; i < 3; i++)
		{
			const auto& pt = cellIter->vertex(i)->point();
			if (bodies[0].frac.has_on_unbounded_side(pt))
				break;

			vert_bounded++;
		}
		if(vert_bounded == 3)
			constrainedCells.push_back(cell_idx);

		cell_idx++;
	}
	/*bool isConstrained;
	for (auto cellIter = cdt.finite_faces_begin(); cellIter != cdt.finite_faces_end(); ++cellIter)
	{
		isConstrained = false;
		for (int i = 0; i < 3; i++)
		{
			isConstrained = cellIter->is_constrained(i);
			if (isConstrained == true)
			{
				constrainedCells.push_back(std::make_pair(cell_idx, i));
				break;
			}
		}
		cell_idx++;
	}*/

	// copy from CDT to IntermediateTriangulation adding info
	// about containing body to triangulation cells
	typedef CDT::Vertex CdtVertex;
	typedef CDT::Face   CdtCell;
	typedef IntermediateTriangulation::Vertex IntermediateVertex;
	typedef IntermediateTriangulation::Face   IntermediateCell;

	DefaultVertexConverter<CdtVertex, IntermediateVertex> vertexConverter;
	BodiesCellConverter<CdtCell, IntermediateCell, CDT> cellConverter(&bodies, &cdt);

	IntermediateTriangulation result;
	copyTriangulation(cdt, result, cellConverter, vertexConverter);
	return result;
}
Polygon Cgal2DMesher::makePolygon(const TaskBody::Border& points) {
	/// convert point set to simple CGAL polygon
	assert(points.size() >= 3); // polygon is a closed line

	Polygon polygon;
	for (const auto& p : points)
		polygon.push_back(CgalPoint2(p[0], p[1]));

	assert(polygon.is_simple()); // has not intersections

	return polygon;
}
CgalPoint2 Cgal2DMesher::findInnerPoint(const Polygon& polygon) {
	/// find inner point in polygon
	CgalVector2 alongBorder(polygon.vertex(0), polygon.vertex(1));
	CgalVector2 crossBorder = alongBorder.perpendicular(CGAL::Orientation::CLOCKWISE);
	CgalPoint2 middle = CGAL::midpoint(polygon.vertex(0), polygon.vertex(1));

	CgalPoint2 innerPoint = middle + crossBorder;
	int n = 1; // iteration number
	while (!polygon.has_on_bounded_side(innerPoint)) {
		// watching on the both sides of the edge, getting closer on each iteration
		innerPoint = middle + crossBorder / pow(-2, n);
		assert(n < 20);
		n++;
	}

	return innerPoint;
}