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

auto load_and_init_mash(){
    vtkSmartPointer<vtkPolyData> mesh = load_mesh(mesh_in_path, false);
    attach_center_to_cells(mesh);
    attach_area(mesh);
    return mesh;
}

auto write_mesh(vtkSmartPointer<vtkPolyData> mesh){
    auto writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer->SetFileName(mesh_out_path.c_str());
    writer->SetInputData(mesh);
    writer->Write();
}
void add_grads(vtkSmartPointer<vtkPolyData> mesh){
    attach_f_true_grad(mesh, grad_f_1);
    grad_calculator grad_calc(f_1, 0.25, kernel_2);
    grad_calc.attach_grad(mesh);
}
int main(int argc, char *argv[]) {
    vtkSmartPointer<vtkPolyData> mesh = load_and_init_mash();
    add_grads(mesh);
    write_mesh(mesh);
}
