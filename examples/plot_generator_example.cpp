#include <cstdlib>
#include <iostream>

#include "../include/functions.h"
#include "../include/geometry.h"
#include "../include/kernels.h"
#include "../include/plot_generator.h"
#include "../include/utils.h"
#include "../include/sphere_generator.h"

const std::string mesh_in_path = "meshes/sphere_642.vtk";
const std::string mesh_out_path = "meshes/sphere_642_with_errors.vtk";

std::string compare_all_norms(vtkSmartPointer<vtkPolyData> mesh,
                              std::function<double(Vector3D)> f,
                              std::function<Vector3D(Vector3D)> grad_f,
                              std::string function_name,double epsilon_min = 0.05, double epsilon_max=1,
                              int num_points = 40) {
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
  double sphere_radius = 1;
  Vector3D sphere_center = Vector3D(0, 0, 0);
  double mesh_size = 0.1;
  sphere_generator m_generator(sphere_radius, sphere_center, mesh_size);
  vtkSmartPointer<vtkPolyData> mesh = m_generator.generate_mesh();
  add_grads(mesh, f_5, grad_f_5, 0.2, kernel_4);
  PlotGenerator generator;
  compare_all_norms(mesh, f_5,  grad_f_5, "acoustic problem");
  // generator.generate_and_linealize_plot(mesh, f_2, grad_f_2, kernel_4, 0.001, 3.0, 0.1, 1.5, 100, "cos(theta)");
}