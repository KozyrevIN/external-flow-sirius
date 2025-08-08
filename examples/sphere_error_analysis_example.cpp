#include <cstdlib>
#include <iostream>
#include <iomanip>
#include <vector>
#include <sstream>
#include <limits>
#include <fstream>
#include <algorithm>
#include <omp.h>

#include <vtkCellData.h>

#include "../include/sphere_generator.h"
#include "../include/functions.h"
#include "../include/geometry.h"
#include "../include/kernels.h"
#include "../include/utils.h"
#include "../include/visualizer.h"
#include "../include/plot_generator.h"

struct AnalysisResult {
    double mesh_size;
    double h_max;  // Maximum cell diameter
    size_t triangle_count;
    double optimal_epsilon;
    double optimal_error;
    double estimated_triangles;
};

// Function to calculate maximum cell diameter (h)
double calculateMaxCellDiameter(vtkSmartPointer<vtkPolyData> mesh) {
    double max_h = 0.0;
    
    for (vtkIdType cellId = 0; cellId < mesh->GetNumberOfCells(); ++cellId) {
        vtkCell* cell = mesh->GetCell(cellId);
        if (cell->GetNumberOfPoints() == 3) { // Triangle
            // Get the three vertices of the triangle
            double p0[3], p1[3], p2[3];
            mesh->GetPoint(cell->GetPointId(0), p0);
            mesh->GetPoint(cell->GetPointId(1), p1);
            mesh->GetPoint(cell->GetPointId(2), p2);
            
            // Calculate edge lengths
            double d01 = sqrt(pow(p1[0]-p0[0], 2) + pow(p1[1]-p0[1], 2) + pow(p1[2]-p0[2], 2));
            double d12 = sqrt(pow(p2[0]-p1[0], 2) + pow(p2[1]-p1[1], 2) + pow(p2[2]-p1[2], 2));
            double d20 = sqrt(pow(p0[0]-p2[0], 2) + pow(p0[1]-p2[1], 2) + pow(p0[2]-p2[2], 2));
            
            // Maximum edge length is the diameter for this triangle
            double cell_h = std::max({d01, d12, d20});
            max_h = std::max(max_h, cell_h);
        }
    }
    
    return max_h;
}

// Function to find optimal epsilon for given mesh
std::pair<double, double> findOptimalEpsilon(vtkSmartPointer<vtkPolyData> mesh, 
                                           double epsilon_min, double epsilon_max, int num_points) {
    double best_epsilon = epsilon_min;
    double best_error = std::numeric_limits<double>::max();
    
    for (int i = 0; i < num_points; ++i) {
        double epsilon = epsilon_min + (epsilon_max - epsilon_min) * i / (num_points - 1);
        
        // Create a copy of the mesh for this epsilon test
        vtkSmartPointer<vtkPolyData> test_mesh = vtkSmartPointer<vtkPolyData>::New();
        test_mesh->DeepCopy(mesh);
        
        // Compute gradients with this epsilon
        add_grads(test_mesh, f_2, grad_f_2, epsilon, kernel_2);
        
        // Compute error
        GradientVisualizer test_visualizer(test_mesh);
        double error = test_visualizer.computeErrorNorm(NormType::L2);
        
        if (error < best_error) {
            best_error = error;
            best_epsilon = epsilon;
        }
    }
    
    return {best_epsilon, best_error};
}

int main(int argc, char *argv[]) {
    try {
        std::cout << "=== Sphere Error Analysis with Multiple Mesh Sizes ===" << std::endl;
        
        // Configuration
        const double sphere_radius = 1.0;
        const Vector3D sphere_center(0.0, 0.0, 0.0);
        const double epsilon_min = 0.1;
        const double epsilon_max = 0.5;
        const int epsilon_samples = 40;
        
        // Different mesh sizes for analysis

        std::vector<double> mesh_sizes;
        for (double i = 0.05; i <= 0.5; i+=0.01 ) {
            mesh_sizes.push_back(i);
        }
        std::vector<AnalysisResult> results;
        
        std::cout << "\nGenerating spheres and finding optimal epsilon for " << mesh_sizes.size() 
                  << " different mesh sizes..." << std::endl;
        std::cout << "Configuration: radius=" << sphere_radius 
                  << ", epsilon_range=[" << epsilon_min << ", " << epsilon_max 
                  << "], samples=" << epsilon_samples << ", function=f_2 (cos(θ))" << std::endl;
        
        #pragma omp parallel for schedule(dynamic, 1) num_threads(8)
        for (size_t i = 0; i < mesh_sizes.size(); ++i) {
            double mesh_size = mesh_sizes[i];
            
            std::cout << "Processing mesh " << (i+1) << "/" << mesh_sizes.size() 
                      << " (size=" << mesh_size << ")..." << std::flush;
            
            // Generate and prepare sphere
            sphere_generator generator(sphere_radius, sphere_center, mesh_size);
            size_t estimated_triangles = generator.estimate_triangle_count();
            
            vtkSmartPointer<vtkPolyData> mesh = generator.generate_mesh();
            size_t actual_triangles = mesh->GetNumberOfCells();
            
            attach_center_to_cells(mesh);
            attach_area(mesh);
            
            // Calculate maximum cell diameter (h)
            double h_max = calculateMaxCellDiameter(mesh);
            
            // Apply function values and true gradients
            attach_f(mesh, f_2);
            attach_f_true_grad(mesh, grad_f_2);
            
            // Find optimal epsilon
            auto [optimal_epsilon, optimal_error] = findOptimalEpsilon(mesh, epsilon_min, epsilon_max, epsilon_samples);
            
            // Generate final mesh with optimal epsilon and save visualization
            add_grads(mesh, f_2, grad_f_2, optimal_epsilon, kernel_2);
            
            GradientVisualizer visualizer(mesh);
            std::ostringstream filename;
            filename << "meshes/sphere_analysis_h" << std::fixed << std::setprecision(3) << mesh_size << ".vtp";
            
            visualizer.add_errors_in_mesh("GradError", "ErrorMagnitude");
            //visualizer.saveVisualization(filename.str());
            
            std::cout << " ✓ (" << actual_triangles << " triangles, h_max=" 
                      << std::fixed << std::setprecision(3) << h_max 
                      << ", optimal_ε=" << std::setprecision(3) << optimal_epsilon << ")" << std::endl;
            
            // Store results
            results.push_back({
                mesh_size, 
                h_max,
                actual_triangles, 
                optimal_epsilon,
                optimal_error,
                static_cast<double>(estimated_triangles)
            });
        }
        
        
        // Write results to CSV file
        std::cout << "\nSaving results..." << std::flush;
        
        // Ensure out directory exists
        std::system("mkdir -p out");
        
        std::string csv_filename = "out/epsilon_vs_h_analysis.csv";
        std::ofstream csv_file(csv_filename);
        
        if (csv_file.is_open()) {
            // Write header
            csv_file << "mesh_size,h_max,triangle_count,optimal_epsilon,optimal_error,epsilon_h_ratio\n";
            
            // Write data
            for (const auto& result : results) {
                double epsilon_h_ratio = result.optimal_epsilon / result.h_max;
                csv_file << std::fixed << std::setprecision(6)
                        << result.mesh_size << ","
                        << result.h_max << ","
                        << result.triangle_count << ","
                        << result.optimal_epsilon << ","
                        << std::scientific << std::setprecision(6) << result.optimal_error << ","
                        << std::fixed << std::setprecision(4) << epsilon_h_ratio << "\n";
            }
            
            csv_file.close();
            std::cout << " ✓" << std::endl;
        } else {
            std::cout << " ✗ Could not create CSV file" << std::endl;
        }
        
        
        std::cout << "\nAnalysis complete!" << std::endl;
        std::cout << "Results saved to: " << csv_filename << std::endl;
        std::cout << "Visualizations saved to: meshes/sphere_analysis_h*.vtp" << std::endl;
        
        if (!results.empty()) {
            double min_ratio = results[0].optimal_epsilon / results[0].h_max;
            double max_ratio = min_ratio;
            for (const auto& result : results) {
                double ratio = result.optimal_epsilon / result.h_max;
                min_ratio = std::min(min_ratio, ratio);
                max_ratio = std::max(max_ratio, ratio);
            }
            std::cout << "Optimal ε/h ratio range: " << std::fixed << std::setprecision(2) 
                      << min_ratio << " - " << max_ratio << std::endl;
        }
        
    } catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return EXIT_FAILURE;
    }
    
    return EXIT_SUCCESS;
}