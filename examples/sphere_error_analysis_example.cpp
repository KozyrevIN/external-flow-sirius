#include <algorithm>
#include <atomic> // For thread-safe counters
#include <chrono> // For high-resolution timing
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <omp.h>
#include <sstream>
#include <vector>
#include <numeric>

#include <vtkCellData.h>

#include "../include/functions.h"
#include "../include/geometry.h"
#include "../include/kernels.h"
#include "../include/plot_generator.h"
#include "../include/sphere_generator.h"
#include "../include/utils.h"
#include "../include/visualizer.h"

// Configuration: Function to test
// Available functions: f_1/grad_f_1, f_2/grad_f_2, f_3/grad_f_3, f_4/grad_f_4, f_5/grad_f_5
const auto TEST_FUNCTION = f_5;
const auto TEST_FUNCTION_GRAD = grad_f_5;
const std::string TEST_FUNCTION_NAME = "f_5_acoustic";

// Example: To switch to f_2 (cos(theta)), change the above to:
// const auto TEST_FUNCTION = f_2;
// const auto TEST_FUNCTION_GRAD = grad_f_2;  
// const std::string TEST_FUNCTION_NAME = "f_2_cos_theta";

struct AnalysisResult {
    double mesh_size;
    double h_max; // Maximum cell diameter
    size_t triangle_count;
    double optimal_epsilon;
    double optimal_error;
    double estimated_triangles;
};

// Function to calculate maximum cell diameter (h)
double calculateMaxCellDiameter(vtkSmartPointer<vtkPolyData> mesh) {
    double max_h = 0.0;

    for (vtkIdType cellId = 0; cellId < mesh->GetNumberOfCells(); ++cellId) {
        vtkCell *cell = mesh->GetCell(cellId);
        if (cell->GetNumberOfPoints() == 3) { // Triangle
            // Get the three vertices of the triangle
            double p0[3], p1[3], p2[3];
            mesh->GetPoint(cell->GetPointId(0), p0);
            mesh->GetPoint(cell->GetPointId(1), p1);
            mesh->GetPoint(cell->GetPointId(2), p2);

            // Calculate edge lengths
            double d01 = sqrt(pow(p1[0] - p0[0], 2) + pow(p1[1] - p0[1], 2) +
                              pow(p1[2] - p0[2], 2));
            double d12 = sqrt(pow(p2[0] - p1[0], 2) + pow(p2[1] - p1[1], 2) +
                              pow(p2[2] - p1[2], 2));
            double d20 = sqrt(pow(p0[0] - p2[0], 2) + pow(p0[1] - p2[1], 2) +
                              pow(p0[2] - p2[2], 2));

            // Maximum edge length is the diameter for this triangle
            double cell_h = std::max({d01, d12, d20});
            max_h = std::max(max_h, cell_h);
        }
    }

    return max_h;
}

// Function to find optimal epsilon for given mesh
std::pair<double, double> findOptimalEpsilon(vtkSmartPointer<vtkPolyData> mesh,
                                             double epsilon_min,
                                             double epsilon_max,
                                             int num_points) {
    double best_epsilon = epsilon_min;
    double best_error = std::numeric_limits<double>::max();

    for (int i = 0; i < num_points; ++i) {
        double epsilon =
            epsilon_min + (epsilon_max - epsilon_min) * i / (num_points - 1);

        // Create a copy of the mesh for this epsilon test
        vtkSmartPointer<vtkPolyData> test_mesh =
            vtkSmartPointer<vtkPolyData>::New();
        test_mesh->DeepCopy(mesh);

        // Compute gradients with this epsilon
        add_grads(test_mesh, TEST_FUNCTION, TEST_FUNCTION_GRAD, epsilon, kernel_2);

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
        std::cout << "=== Sphere Error Analysis with Multiple Mesh Sizes ==="
                  << std::endl;

        // Configuration
        const double sphere_radius = 1.0;
        const Vector3D sphere_center(0.0, 0.0, 0.0);
        const double epsilon_min = 0.1;
        const double epsilon_max = 0.5;
        const int epsilon_samples = 40;

        // Different mesh sizes for analysis  
        std::vector<double> mesh_sizes;
        for (double i = 0.02; i <= 0.25; i += 0.01) {  // Original range for full analysis
            mesh_sizes.push_back(i);
        }
        std::vector<AnalysisResult> results;
        results.reserve(mesh_sizes.size()); // Pre-allocate memory

        // --- Timing variables ---
        double time_mesh_generation = 0.0;
        double time_data_attachment = 0.0;
        double time_epsilon_search = 0.0;
        double time_other_and_storage = 0.0;

        std::cout << "\nProcessing " << mesh_sizes.size() << " different mesh sizes using up to 8 threads..." << std::flush;
        auto loop_start_time = std::chrono::high_resolution_clock::now();

        // The reduction clause creates a private copy of each timer for each thread,
        // and safely sums them all together at the end of the loop.
#pragma omp parallel for schedule(dynamic, 1) num_threads(8) \
    reduction(+:time_mesh_generation, time_data_attachment, time_epsilon_search, time_other_and_storage)
        for (size_t i = 0; i < mesh_sizes.size(); ++i) {
            auto last_timestamp = std::chrono::high_resolution_clock::now();
            
            double mesh_size = mesh_sizes[i];

            // --- Section 1: Mesh Generation ---
            sphere_generator generator(sphere_radius, sphere_center, mesh_size);
            vtkSmartPointer<vtkPolyData> mesh;
#pragma omp critical(mesh_generation) // Named critical section for clarity
            {
                mesh = generator.generate_mesh();
            }
            auto timestamp_after_mesh = std::chrono::high_resolution_clock::now();
            time_mesh_generation += std::chrono::duration<double>(timestamp_after_mesh - last_timestamp).count();
            last_timestamp = timestamp_after_mesh;

            // --- Section 2: Data Attachment and Mesh Prep ---
            size_t actual_triangles = mesh->GetNumberOfCells();
            attach_center_to_cells(mesh);
            attach_area(mesh);
            double h_max = calculateMaxCellDiameter(mesh);
            attach_f(mesh, TEST_FUNCTION);
            attach_f_true_grad(mesh, TEST_FUNCTION_GRAD);
            auto timestamp_after_attach = std::chrono::high_resolution_clock::now();
            time_data_attachment += std::chrono::duration<double>(timestamp_after_attach - last_timestamp).count();
            last_timestamp = timestamp_after_attach;
            
            // --- Section 3: Find Optimal Epsilon ---
            auto [optimal_epsilon, optimal_error] = findOptimalEpsilon(
                mesh, epsilon_min, epsilon_max, epsilon_samples);
            auto timestamp_after_epsilon = std::chrono::high_resolution_clock::now();
            time_epsilon_search += std::chrono::duration<double>(timestamp_after_epsilon - last_timestamp).count();
            last_timestamp = timestamp_after_epsilon;
            
            // --- Section 3.5: Compute final gradients with optimal epsilon and save mesh ---
#pragma omp critical(mesh_processing_and_saving) // Thread-safe processing and file writing
            {
                // Create final mesh with optimal epsilon
                vtkSmartPointer<vtkPolyData> final_mesh = vtkSmartPointer<vtkPolyData>::New();
                final_mesh->DeepCopy(mesh);
                add_grads(final_mesh, TEST_FUNCTION, TEST_FUNCTION_GRAD, optimal_epsilon, kernel_2);
                
                // Add error visualization to the mesh
                GradientVisualizer visualizer(final_mesh);
                visualizer.add_errors_in_mesh("GradientDifference", "ErrorMagnitude");
                
                // Save mesh with all fields attached
                std::ostringstream mesh_filename;
                mesh_filename << "meshes/sphere_analysis_" << TEST_FUNCTION_NAME 
                             << "_h" << std::fixed << std::setprecision(3) 
                             << mesh_size << "_eps" << std::setprecision(3) << optimal_epsilon << ".vtp";
                write_mesh(final_mesh, mesh_filename.str());
                std::cout << "Saved mesh: " << mesh_filename.str() << std::endl;
            }

            // --- Section 4: Storage (and any other work) ---
#pragma omp critical(results_storage) // Named critical section
            {
                results.push_back({mesh_size, h_max, actual_triangles,
                                   optimal_epsilon, optimal_error,
                                   static_cast<double>(generator.estimate_triangle_count())});
            }
            auto timestamp_after_storage = std::chrono::high_resolution_clock::now();
            time_other_and_storage += std::chrono::duration<double>(timestamp_after_storage - last_timestamp).count();
        }

        auto loop_end_time = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> total_loop_duration = loop_end_time - loop_start_time;
        std::cout << " ✓ Done in " << std::fixed << std::setprecision(2) << total_loop_duration.count() << " seconds." << std::endl;

        // --- Timing Report ---
        double total_cpu_time = time_mesh_generation + time_data_attachment + time_epsilon_search + time_other_and_storage;
        std::cout << "\n--- Performance Analysis ---\n";
        std::cout << std::fixed << std::setprecision(4);
        std::cout << "Total Wall-Clock Time: " << total_loop_duration.count() << " s\n";
        std::cout << "Total CPU Time (sum of all threads): " << total_cpu_time << " s\n\n";
        std::cout << "Breakdown of CPU time:\n";
        std::cout << "  - Mesh Generation:    " << std::setw(8) << time_mesh_generation << " s (" << std::setw(5) << std::setprecision(1) << (time_mesh_generation / total_cpu_time * 100.0) << "%)\n";
        std::cout << "  - Data Attachment:    " << std::setw(8) << time_data_attachment << " s (" << std::setw(5) << std::setprecision(1) << (time_data_attachment / total_cpu_time * 100.0) << "%)\n";
        std::cout << "  - Epsilon Search:     " << std::setw(8) << time_epsilon_search << " s (" << std::setw(5) << std::setprecision(1) << (time_epsilon_search / total_cpu_time * 100.0) << "%)\n";
        std::cout << "  - Storage & Other:    " << std::setw(8) << time_other_and_storage << " s (" << std::setw(5) << std::setprecision(1) << (time_other_and_storage / total_cpu_time * 100.0) << "%)\n";
        std::cout << "---------------------------\n";


        // Write results to CSV file
        std::cout << "\nSaving results..." << std::flush;
        std::system("mkdir -p out");
        std::string csv_filename = "out/epsilon_vs_h_analysis_" + TEST_FUNCTION_NAME + ".csv";
        std::ofstream csv_file(csv_filename);

        if (csv_file.is_open()) {
            csv_file << "mesh_size,h_max,triangle_count,optimal_epsilon,optimal_error,epsilon_h_ratio\n";
            for (const auto &result : results) {
                double epsilon_h_ratio = result.optimal_epsilon / result.h_max;
                csv_file << std::fixed << std::setprecision(6)
                         << result.mesh_size << "," << result.h_max << ","
                         << result.triangle_count << ","
                         << result.optimal_epsilon << "," << std::scientific
                         << std::setprecision(6) << result.optimal_error << ","
                         << std::fixed << std::setprecision(4)
                         << epsilon_h_ratio << "\n";
            }
            csv_file.close();
            std::cout << " ✓" << std::endl;
        } else {
            std::cout << " ✗ Could not create CSV file" << std::endl;
        }

        std::cout << "\nAnalysis complete!" << std::endl;
        std::cout << "Results saved to: " << csv_filename << std::endl;

        if (!results.empty()) {
            double min_ratio = results[0].optimal_epsilon / results[0].h_max;
            double max_ratio = min_ratio;
            for (const auto &result : results) {
                double ratio = result.optimal_epsilon / result.h_max;
                min_ratio = std::min(min_ratio, ratio);
                max_ratio = std::max(max_ratio, ratio);
            }
            std::cout << "Optimal ε/h ratio range: " << std::fixed
                      << std::setprecision(2) << min_ratio << " - " << max_ratio
                      << std::endl;
        }

    } catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}