#include <cstdlib>
#include <iostream>
#include <string>


#include <vtkCellData.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkXMLPolyDataWriter.h>


#include "../include/functions.h"
#include "../include/geometry.h"
#include "../include/kernels.h"
#include "../include/utils.h"
#include "../include/visualizer.h"


const std::string mesh_in_path = "meshes/sphere_642.vtk";
const std::string mesh_out_path = "meshes/sphere_642_w_fields.vtp";

int main(int argc, char *argv[]) {
    vtkSmartPointer<vtkPolyData> mesh = load_and_init_mash(mesh_in_path);
    add_grads(mesh, f_2, grad_f_2, 0.25, kernel_2);
    write_mesh(mesh, mesh_out_path);
}
