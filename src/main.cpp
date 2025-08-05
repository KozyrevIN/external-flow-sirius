#include <iostream>
#include <string>

#include <vtkPolyData.h>
#include <vtkSmartPointer.h>

#include "../include/utils.h"
#include "../include/kernels.h"
const std::string mesh_path = "meshes/sphere_642.vtk";
int main(int argc, char *argv[]) {
    std::cout << kernel_2(1.0) << std::endl;
    return EXIT_SUCCESS;
}
