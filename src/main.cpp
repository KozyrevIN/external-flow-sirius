#include <iostream>
#include <string>

#include <vtkPolyData.h>
#include <vtkSmartPointer.h>

#include "../include/utils.h"
const std::string mesh_path = "meshes/sphere_642.vtk";
int main(int argc, char *argv[]) {
    vtkSmartPointer<vtkPolyData> mesh = load_mesh(mesh_path, true);

    return EXIT_SUCCESS;
}
