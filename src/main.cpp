#include <iostream>
#include <string>

#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkXMLPolyDataWriter.h>

#include "../include/functions.h"
#include "../include/kernels.h"
#include "../include/utils.h"
#include "../include/geometry.h"

const std::string mesh_in_path = "meshes/sphere_642.vtk";
const std::string mesh_out_path = "meshes/sphere_642_w_fields.vtp";

int main(int argc, char *argv[]) {
    auto mesh = load_mesh(mesh_in_path, false);

    add_center_to_cells(mesh);
    attach_f(mesh, f_2);
    attach_area(mesh);
    attach_f_true_grad(mesh, grad_f_2);

    auto writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer->SetFileName(mesh_out_path.c_str());
    writer->SetInputData(mesh);
    writer->Write();

    return EXIT_SUCCESS;
}
