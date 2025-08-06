#include <iostream>
#include <string>

#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkXMLPolyDataWriter.h>

#include "../include/attach_f.h"
#include "../include/functions.h"
#include "../include/kernels.h"
#include "../include/utils.h"

const std::string mesh_in_path = "meshes/sphere_642.vtk";
const std::string mesh_out_path = "meshes/sphere_642_w_fields.vtp";

int main(int argc, char *argv[]) {
    auto mesh = load_mesh(mesh_in_path, false);

    attach_f(mesh, f_2);

    auto writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer->SetFileName(mesh_out_path.c_str());
    writer->SetInputData(mesh);
    writer->Write();

    return EXIT_SUCCESS;
}
