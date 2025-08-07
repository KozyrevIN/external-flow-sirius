#include <iostream>
#include <string>
#include <cstdlib>

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
    vtkSmartPointer<vtkPolyData> mesh = load_mesh(mesh_in_path, false);

    // Сначала создаем необходимые атрибуты
    attach_center_to_cells(mesh);
    attach_area(mesh);
    compute_f_true_grad(mesh, grad_f_2);
    attach_f(mesh, f_2);

    grad_calculator grad_calc(f_2, 0.1, kernel_2);
    grad_calc.attach_grad(mesh);

    auto writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer->SetFileName(mesh_out_path.c_str());
    writer->SetInputData(mesh);
    writer->Write();
    return EXIT_SUCCESS;
}
