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

#include "VTKSnapshotter.hpp"

#include "Model.hpp"

using std::vector;

template <class modelType>
const std::string VTKSnapshotter<modelType>::prefix = "";

template<class modelType>
VTKSnapshotter<modelType>::VTKSnapshotter(const Mesh* _mesh) : mesh(_mesh)
{
	pattern = prefix + "CGA_First_%{STEP}.vtu";
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
void VTKSnapshotter<FirstModel>::dump(const int i)
{
	auto grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
	auto points = vtkSmartPointer<vtkPoints>::New();

	for (const auto& vhandle : mesh->vertexHandles)
	{
		const auto& pt = vhandle->point();
		points->InsertNextPoint(pt[0], pt[1], 0.0);
	}
	grid->SetPoints(points);

	auto tris= vtkSmartPointer<vtkCellArray>::New();
	for (const auto& cell : mesh->cellHandles) 
	{
		auto vtkCell = vtkSmartPointer<vtkTriangle>::New();

		for (int i = 0; i < Mesh::CELL_POINTS_NUMBER; i++) 
			vtkCell->GetPointIds()->SetId(i, cell->info().localVertexIndices[i]);

		tris->InsertNextCell(vtkCell);
	}

	grid->SetCells(VTK_TRIANGLE, tris);
	vtkCellData* fd = grid->GetCellData();

	// Writing
	auto writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
	writer->SetFileName(getFileName(i).c_str());
	writer->SetInputData(grid);
	writer->Write();
}

template class VTKSnapshotter<FirstModel>;