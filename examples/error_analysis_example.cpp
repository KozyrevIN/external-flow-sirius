#include <cstdlib>
#include <iostream>

#include <vtkCellData.h>

#include "../include/functions.h"
#include "../include/geometry.h"
#include "../include/kernels.h"
#include "../include/utils.h"
#include "../include/visualizer.h"
#include "../include/plot_generator.h"

const std::string mesh_in_path = "meshes/sphere_642.vtk";

int main(int argc, char *argv[]) {
    vtkSmartPointer<vtkPolyData> mesh = load_and_init_mash(mesh_in_path);
    
    add_grads(mesh, f_3, grad_f_3, 0.05, kernel_2);
    
    GradientVisualizer visualizer(mesh);
    
    visualizer.add_errors_in_mesh("GradError", "ErrorMagnitude");
    
    visualizer.saveVisualization("error_analysis.vtp");   
}