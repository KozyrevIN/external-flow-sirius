#include "../include/utils.h"

#include <iostream>
#include <vtkCellData.h>
#include <vtkDataSetSurfaceFilter.h>
#include <vtkGeometryFilter.h>
#include <vtkPolyDataWriter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkXMLPolyDataWriter.h>

#include "../include/functions.h"
#include "../include/geometry.h"


vtkSmartPointer<vtkPolyData> load_mesh(std::string filename, bool verbose) {

    // Create the reader for UnstructuredGrid data
    vtkSmartPointer<vtkUnstructuredGridReader> reader =
        vtkSmartPointer<vtkUnstructuredGridReader>::New();
    reader->SetFileName(filename.c_str());
    reader->Update();

    // Get the output and cast it to an UnstructuredGrid
    vtkSmartPointer<vtkUnstructuredGrid> grid =
        vtkUnstructuredGrid::SafeDownCast(reader->GetOutput());

    if (verbose) {
        std::cout << "Successfully loaded " << filename << std::endl;
        std::cout << "------------------------------------------" << std::endl;
        std::cout << "Number of points (vertices): "
                  << grid->GetNumberOfPoints() << std::endl;
        std::cout << "Number of cells (triangles): " << grid->GetNumberOfCells()
                  << std::endl;
        std::cout << "------------------------------------------" << std::endl;
    }

    // Use DataSetSurfaceFilter instead of GeometryFilter to avoid memory leak
    auto surfaceFilter = vtkSmartPointer<vtkDataSetSurfaceFilter>::New();
    surfaceFilter->SetInputData(grid);
    surfaceFilter->Update();

    // Create a deep copy of the output
    vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
    polydata->DeepCopy(surfaceFilter->GetOutput());

    // Clean up filter connections
    surfaceFilter->SetInputData(nullptr);

    return polydata;
}

void save_mesh(vtkPolyData *polyData, std::string filepath) {
    vtkSmartPointer<vtkPolyDataWriter> writer =
        vtkSmartPointer<vtkPolyDataWriter>::New();
    writer->SetFileName(filepath.c_str());
    writer->SetInputData(polyData);
    writer->Write();

    std::cout << "Mesh saved to " << filepath << std::endl;
}

Vector3D getCenter(vtkIdType cellId, vtkPolyData *polyData) {
    vtkDataArray *centerArray =
        polyData->GetCellData()->GetArray("CellCenters");

    double center[3];
    centerArray->GetTuple(cellId, center);
    Vector3D center_vector = Vector3D(center[0], center[1], center[2]);
    return center_vector;
}

double getAttributeArea(vtkIdType cellId, vtkPolyData *polyData) {
    vtkDataArray *areaArray = polyData->GetCellData()->GetArray("Area");
    return areaArray->GetTuple1(cellId);
}

vtkSmartPointer<vtkPolyData> load_and_init_mash(std::string mesh_in_path){
    vtkSmartPointer<vtkPolyData> mesh = load_mesh(mesh_in_path, false);
    attach_center_to_cells(mesh);
    attach_area(mesh);
    return mesh;
}

void write_mesh(vtkSmartPointer<vtkPolyData> mesh, std::string mesh_out_path){
    // Ensure meshes directory exists
    std::string mkdir_cmd = "mkdir -p meshes";
    std::system(mkdir_cmd.c_str());
    
    // Prepend meshes/ if not already present
    if (mesh_out_path.find("meshes/") != 0) {
        mesh_out_path = "meshes/" + mesh_out_path;
    }
    
    auto writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer->SetFileName(mesh_out_path.c_str());
    writer->SetInputData(mesh);
    writer->Write();
}
void add_grads(vtkSmartPointer<vtkPolyData> mesh, std::function<double(Vector3D)> f, std::function<Vector3D(Vector3D)> grad_f, double epsilon, std::function<double(double)> kernel){
    attach_f_true_grad(mesh, grad_f);
    grad_calculator grad_calc(f, epsilon, kernel);
    grad_calc.attach_grad(mesh);
}