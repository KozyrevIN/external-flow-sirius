#include <cstdlib>
#include <iostream>

#include "../include/functions.h"
#include "../include/geometry.h"
#include "../include/kernels.h"
#include "../include/plot_generator.h"
#include "../include/utils.h"

const std::string mesh_in_path = "meshes/sphere_642.vtk";
const std::string mesh_out_path = "meshes/sphere_642_with_errors.vtk";


int main(int argc, char *argv[]) {
    vtkSmartPointer<vtkPolyData> mesh = load_and_init_mash(mesh_in_path);
    
    add_grads(mesh, f_2, grad_f_2, 0.25, kernel_2);
    
    PlotGenerator generator;
    
    std::string comparison_plot = generator.generateKernelComparisonPlot(
        mesh, f_2, grad_f_2, 0.001, 10.0, 100, "cos(theta)");
    
    std::cout << "Generated comparison plot: " << comparison_plot << std::endl;

}