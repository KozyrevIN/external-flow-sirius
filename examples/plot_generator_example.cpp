#include <cstdlib>
#include <iostream>

#include "../include/functions.h"
#include "../include/geometry.h"
#include "../include/kernels.h"
#include "../include/plot_generator.h"
#include "../include/utils.h"

const std::string mesh_in_path = "meshes/sphere_642.vtk";
const std::string mesh_out_path = "meshes/sphere_642_with_errors.vtk";

std::string compare_all_norms(vtkSmartPointer<vtkPolyData> mesh,
                              std::function<double(Vector3D)> f,
                              std::function<Vector3D(Vector3D)> grad_f,
                              std::string function_name,double epsilon_min = 0.01, double epsilon_max=5,
                              int num_points = 100) {
  PlotGenerator generator;
  std::string comparison_plot = generator.generateKernelComparisonPlot(
      mesh, f, grad_f, epsilon_min, epsilon_max, num_points, function_name, NormType::LINF);
  std::string comparison_plot2 = generator.generateKernelComparisonPlot(
      mesh, f, grad_f, epsilon_min, epsilon_max, num_points, function_name, NormType::L2);
  std::string comparison_plot3 = generator.generateKernelComparisonPlot(
      mesh, f, grad_f, epsilon_min, epsilon_max, num_points, function_name, NormType::L1);
  return comparison_plot + comparison_plot2 + comparison_plot3;
}

int main(int argc, char *argv[]) {
  vtkSmartPointer<vtkPolyData> mesh = load_and_init_mash(mesh_in_path);
  add_grads(mesh, f_2, grad_f_2, 0.25, kernel_2);
  PlotGenerator generator;
  
  // compare_all_norms(mesh, f_4,  grad_f_4, "exp(y*y-z*z)");
  generator.generateKernelComparisonPlot(mesh, f_2, grad_f_2, 0.001, 3.0, 100, "cos(theta)", NormType::L2);
  // generator.generate_and_linealize_plot(mesh, f_4, grad_f_4, kernel_4, 0.1, 2, 0.15, 1, 100, "exp(y*y-z*z)");
}