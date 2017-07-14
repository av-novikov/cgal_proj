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
#include "src/models/Acid/Acid2d.hpp"

using std::vector;

template<class modelType>
const std::string VTKSnapshotter<modelType>::prefix = "snaps/";
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
	auto id = vtkSmartPointer<vtkIntArray>::New();
	id->SetName("id");

	auto pres = vtkSmartPointer<vtkDoubleArray>::New();
	pres->SetName("pressure");
	auto perm_x = vtkSmartPointer<vtkDoubleArray>::New();
	perm_x->SetName("k_x");
	auto perm_y = vtkSmartPointer<vtkDoubleArray>::New();
	perm_y->SetName("k_y");

	points->Allocate(mesh->getVerticesSize());
	facets->Allocate(mesh->getCellsSize());

	for (const auto& vhandle : mesh->vertexHandles)
	{
		const auto& pt = vhandle->point();
		points->InsertNextPoint(pt[0] * model->R_dim, pt[1] * model->R_dim, 0.0);
	}
	for (int i = 0; i < mesh->inner_cells; i++) 
	{
		const Cell& cell = mesh->cells[i];
		auto vtkCell = vtkSmartPointer<vtkTriangle>::New();

		for (int i = 0; i < Mesh::CELL_POINTS_NUMBER; i++) 
			vtkCell->GetPointIds()->SetId(i, cell.points[i]);

		facets->InsertNextCell(vtkCell);
		type->InsertNextValue(cell.type);
		vol->InsertNextValue(cell.V);
		id->InsertNextValue(cell.id);
		const auto& var = (*model)[i].u_next;
		pres->InsertNextValue(var.p * model->P_dim / BAR_TO_PA);
		perm_x->InsertNextValue(model->getPerm(cell) * model->R_dim * model->R_dim);
		perm_y->InsertNextValue(model->getPerm(cell) * model->R_dim * model->R_dim);
	}

	grid->SetPoints(points);
	grid->SetCells(VTK_TRIANGLE, facets);
	vtkCellData* fd = grid->GetCellData();
	fd->AddArray(type);
	fd->AddArray(vol);
	fd->AddArray(id);
	fd->AddArray(pres);
	fd->AddArray(perm_x);
	fd->AddArray(perm_y);

	// Writing
	auto writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
	writer->SetFileName(getFileName(i).c_str());
	writer->SetInputData(grid);
	writer->Write();
}
void VTKSnapshotter<acid2d::Acid2d>::dump(const int i)
{
	using namespace acid2d;
	auto grid = vtkSmartPointer<vtkUnstructuredGrid>::New();
	auto points = vtkSmartPointer<vtkPoints>::New();
	auto facets = vtkSmartPointer<vtkCellArray>::New();

	auto type = vtkSmartPointer<vtkIntArray>::New();
	type->SetName("type");
	auto poro = vtkSmartPointer<vtkDoubleArray>::New();
	poro->SetName("porosity");
	auto pres = vtkSmartPointer<vtkDoubleArray>::New();
	pres->SetName("pressure");
	auto s_wat = vtkSmartPointer<vtkDoubleArray>::New();
	s_wat->SetName("s_wat");
	auto s_oil = vtkSmartPointer<vtkDoubleArray>::New();
	s_oil->SetName("s_oil");
	auto xa = vtkSmartPointer<vtkDoubleArray>::New();
	xa->SetName("x_acid");
	auto xw = vtkSmartPointer<vtkDoubleArray>::New();
	xw->SetName("x_water");
	auto xs = vtkSmartPointer<vtkDoubleArray>::New();
	xs->SetName("x_salt");
	auto perm_x = vtkSmartPointer<vtkDoubleArray>::New();
	perm_x->SetName("k_x");
	auto perm_y = vtkSmartPointer<vtkDoubleArray>::New();
	perm_y->SetName("k_y");

	points->Allocate(mesh->getVerticesSize());
	facets->Allocate(mesh->getCellsSize());

	for (const auto& vhandle : mesh->vertexHandles)
	{
		const auto& pt = vhandle->point();
		points->InsertNextPoint(pt[0] * model->R_dim, pt[1] * model->R_dim, 0.0);
	}
	for (int i = 0; i < mesh->inner_cells; i++)
	{
		const Cell& cell = mesh->cells[i];
		auto vtkCell = vtkSmartPointer<vtkTriangle>::New();

		for (int i = 0; i < Mesh::CELL_POINTS_NUMBER; i++)
			vtkCell->GetPointIds()->SetId(i, cell.points[i]);

		facets->InsertNextCell(vtkCell);
		type->InsertNextValue(cell.type);
		const auto& var = (*model)[i].u_next;
		poro->InsertNextValue(var.m);
		pres->InsertNextValue(var.p * model->P_dim / BAR_TO_PA);
		s_wat->InsertNextValue(var.s);
		s_oil->InsertNextValue(1.0 - var.s);
		xa->InsertNextValue(var.xa);
		xw->InsertNextValue(var.xw);
		xs->InsertNextValue(1.0 - var.xa - var.xw);
		perm_x->InsertNextValue(M2toMilliDarcy(model->getPermValue(cell) * model->R_dim * model->R_dim));
		perm_y->InsertNextValue(M2toMilliDarcy(model->getPerm(cell).value() * model->R_dim * model->R_dim));
	}

	grid->SetPoints(points);
	grid->SetCells(VTK_TRIANGLE, facets);
	vtkCellData* fd = grid->GetCellData();
	fd->AddArray(type);
	fd->AddArray(poro);
	fd->AddArray(pres);
	fd->AddArray(s_oil);
	fd->AddArray(s_wat);
	fd->AddArray(xa);
	fd->AddArray(xw);
	fd->AddArray(xs);
	fd->AddArray(perm_x);
	fd->AddArray(perm_y);

	// Writing
	auto writer = vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
	writer->SetFileName(getFileName(i).c_str());
	writer->SetInputData(grid);
	writer->Write();
}

template class VTKSnapshotter<oil2d::Oil2d>;
template class VTKSnapshotter<acid2d::Acid2d>;