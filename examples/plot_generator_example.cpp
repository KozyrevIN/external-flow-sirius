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
                              double epsilon_min, double epsilon_max,
                              int num_points, std::string function_name) {
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
  add_grads(mesh, f_3, grad_f_3, 0.25, kernel_2);
  PlotGenerator generator;

  std::string comparison_plot = compare_all_norms(mesh, f_3, grad_f_3, 0.0001, 10.0, 100, "x^2+y^2+z^2");
  std::cout << "Generated comparison plot: " << comparison_plot << std::endl;
}