#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkDoubleArray.h>
#include <vtkIntArray.h>
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkPoints.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkTriangle.h>
#include <vtkTetra.h>

#include "src/snapshotter/VTKSnapshotter.hpp"
#include "src/models/Oil2d/Oil2d.hpp"

using std::vector;

template<class modelType>
const std::string VTKSnapshotter<modelType>::prefix = "";
template<class modelType>
VTKSnapshotter<modelType>::VTKSnapshotter(const modelType* _model) : model(_model), mesh(_model->getMesh())
{
	pattern = prefix + "CGAL_First_%{STEP}.vtu";

}
template<class modelType>
VTKSnapshotter<modelType>::~VTKSnapshotter()
{
}
template<class modelType>
std::string VTKSnapshotter<modelType>::getFileName(int i)
{
	std::string filename = pattern;
	return replace(filename, "%{STEP}", std::to_string(i));
};
template<class modelType>
std::string VTKSnapshotter<modelType>::replace(std::string filename, std::string from, std::string to)
{
	size_t start_pos = 0;
	while ((start_pos = filename.find(from, start_pos)) != std::string::npos)
	{
		filename.replace(start_pos, from.length(), to);
		start_pos += to.length();
	}
	return filename;
};

template<class modelType>
void VTKSnapshotter<modelType>::dump(const int i)
{
}
void VTKSnapshotter<oil2d::Oil2d>::dump(const int i)
{
	using namespace oil2d;
	auto grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
	auto points = vtkSmartPointer<vtkPoints>::New();
	auto facets = vtkSmartPointer<vtkCellArray>::New();
	auto type = vtkSmartPointer<vtkIntArray>::New();
	auto vol = vtkSmartPointer<vtkDoubleArray>::New();
	type->SetName("type");
	vol->SetName("volume");

	points->Allocate(mesh->triangulation.number_of_vertices());
	facets->Allocate(mesh->triangulation.number_of_faces());

	for (const auto& vhandle : mesh->vertexHandles)
	{
		const auto& pt = vhandle->point();
		points->InsertNextPoint(pt[0], pt[1], 0.0);
	}
	for (int i = 0; i < mesh->inner_cells; i++) 
	{
		const TriangleCell& cell = mesh->cells[i];
		auto vtkCell = vtkSmartPointer<vtkTriangle>::New();

		for (int i = 0; i < Mesh::CELL_POINTS_NUMBER; i++) 
			vtkCell->GetPointIds()->SetId(i, cell.points[i]);

		facets->InsertNextCell(vtkCell);
		type->InsertNextValue(cell.type);
		vol->InsertNextValue(cell.V);
	}

	grid->SetPoints(points);
	grid->SetCells(VTK_TRIANGLE, facets);
	vtkCellData* fd = grid->GetCellData();
	fd->AddArray(type);
	fd->AddArray(vol);

	// Writing
	auto writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
	writer->SetFileName(getFileName(i).c_str());
	writer->SetInputData(grid);
	writer->Write();
}

template class VTKSnapshotter<oil2d::Oil2d>;