#include "../include/utils.h"

#include <vtkGeometryFilter.h>
#include <vtkUnstructuredGrid.h>
#include <vtkUnstructuredGridReader.h>
#include <vtkPolyDataWriter.h>
#include <vtkCellData.h>
#include <iostream>

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

    // This filter extracts the surface geometry from any dataset.
    auto geometryFilter = vtkSmartPointer<vtkGeometryFilter>::New();
    geometryFilter->SetInputData(grid);
    geometryFilter->Update();

    // The output of the filter is the vtkPolyData
    vtkSmartPointer<vtkPolyData> polydata = geometryFilter->GetOutput();

    return polydata;
}

void save_mesh(vtkPolyData* polyData, std::string filepath) {
    vtkSmartPointer<vtkPolyDataWriter> writer = vtkSmartPointer<vtkPolyDataWriter>::New();
    writer->SetFileName(filepath.c_str());
    writer->SetInputData(polyData);
    writer->Write();
    
    std::cout << "Mesh saved to " << filepath << std::endl;
}

Vector3D getCenter(vtkIdType cellId, vtkPolyData* polyData) {
    vtkDataArray* centerArray = polyData->GetCellData()->GetArray("CellCenters");
    
    double center[3];
    centerArray->GetTuple(cellId, center);
    Vector3D center_vector = Vector3D(center[0], center[1], center[2]);
    return center_vector;
}

double getAttributeArea(vtkIdType cellId, vtkPolyData* polyData) {
    vtkDataArray* areaArray = polyData->GetCellData()->GetArray("Area");
    return areaArray->GetTuple1(cellId);
}
