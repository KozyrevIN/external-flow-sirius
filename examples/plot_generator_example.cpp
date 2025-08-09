#include <cstdlib>
#include <iostream>
#include <chrono>

#include "../include/functions.h"
#include "../include/geometry.h"
#include "../include/kernels.h"
#include "../include/plot_generator.h"
#include "../include/utils.h"
#include "../include/sphere_generator.h"
#include "../include/visualizer.h"

// Configuration: Function to test and analysis parameters
// Available functions: f_1/grad_f_1, f_2/grad_f_2, f_3/grad_f_3, f_4/grad_f_4, f_5/grad_f_5
const auto TEST_FUNCTION = f_5;
const auto TEST_FUNCTION_GRAD = grad_f_5;
const std::string TEST_FUNCTION_NAME = "f_5_acoustic";

// Analysis parameters
const double TEST_EPSILON = 0.2;
const auto TEST_KERNEL = kernel_4;
const std::string TEST_KERNEL_NAME = "kernel_4";

// Example: To switch to f_2 (cos(theta)), change the above to:
// const auto TEST_FUNCTION = f_2;
// const auto TEST_FUNCTION_GRAD = grad_f_2;  
// const std::string TEST_FUNCTION_NAME = "f_2_cos_theta";

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
  double mesh_size = 0.05;  // Larger mesh size for faster testing of parallelization
  sphere_generator m_generator(sphere_radius, sphere_center, mesh_size);
  vtkSmartPointer<vtkPolyData> mesh = m_generator.generate_mesh();
  
  // Initialize geometric properties (centers and areas) before using gradient functions
  attach_center_to_cells(mesh);
  attach_area(mesh);
  attach_f(mesh, TEST_FUNCTION);
  
  add_grads(mesh, TEST_FUNCTION, TEST_FUNCTION_GRAD, TEST_EPSILON, TEST_KERNEL);
  GradientVisualizer vis(mesh, "f true grad", "Grad");
  vis.add_errors_in_mesh("GradientDifference", "ErrorMagnitude");
  
  // Generate a single comparison plot (testing parallelization)
  std::cout << "Starting kernel comparison with parallel epsilon evaluation..." << std::endl;
  auto start_time = std::chrono::high_resolution_clock::now();
  
  PlotGenerator generator;
  std::string plot_result = generator.generateKernelComparisonPlot(
      mesh, TEST_FUNCTION, TEST_FUNCTION_GRAD, 0.05, 1.0, 30, TEST_FUNCTION_NAME, NormType::L2);
  
  auto end_time = std::chrono::high_resolution_clock::now();
  auto duration = std::chrono::duration_cast<std::chrono::seconds>(end_time - start_time);
  
  std::cout << "Generated kernel comparison plot: " << plot_result << std::endl;
  std::cout << "Parallel kernel comparison completed in " << duration.count() << " seconds" << std::endl;
  
  // Save mesh with all computed fields
  std::string mesh_out_path = "meshes/plot_generator_" + TEST_FUNCTION_NAME + 
                              "_eps" + std::to_string(TEST_EPSILON) + 
                              "_" + TEST_KERNEL_NAME + ".vtp";
  write_mesh(mesh, mesh_out_path);
  
  std::cout << "\n=== RESULTS SUMMARY ===" << std::endl;
  std::cout << "Function tested: " << TEST_FUNCTION_NAME << std::endl;
  std::cout << "Epsilon used: " << TEST_EPSILON << std::endl;
  std::cout << "Kernel used: " << TEST_KERNEL_NAME << std::endl;
  std::cout << "Mesh cells: " << mesh->GetNumberOfCells() << std::endl;
  std::cout << "Mesh saved with fields:" << std::endl;
  std::cout << "  - CellCenters (3D coordinates)" << std::endl;
  std::cout << "  - Area (cell surface areas)" << std::endl;
  std::cout << "  - f (function values)" << std::endl;
  std::cout << "  - f true grad (analytical gradients)" << std::endl;
  std::cout << "  - Grad (computed gradients)" << std::endl;
  std::cout << "  - GradientDifference (error vectors)" << std::endl;
  std::cout << "  - ErrorMagnitude (error magnitudes)" << std::endl;
  std::cout << "Mesh saved to: " << mesh_out_path << std::endl;
}