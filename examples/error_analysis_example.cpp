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
    try {
        // Load and prepare mesh
        vtkSmartPointer<vtkPolyData> mesh = load_and_init_mash(mesh_in_path);
        
        // Add gradients for analysis
        add_grads(mesh, f_2, grad_f_2, 0.3, kernel_2);
        
        std::cout << "=== Complete Error Analysis Example ===" << std::endl;
        
        // Create visualizer for analysis (read-only)
        GradientVisualizer visualizer(mesh);
        
        // 1. Compute error statistics without modifying mesh
        std::cout << "\n1. Computing error statistics..." << std::endl;
        double l1_error = visualizer.computeErrorNorm(NormType::L1);
        double l2_error = visualizer.computeErrorNorm(NormType::L2);
        double linf_error = visualizer.computeErrorNorm(NormType::LINF);
        
        std::cout << "Error Statistics:" << std::endl;
        std::cout << "  L1 norm:  " << l1_error << std::endl;
        std::cout << "  L2 norm:  " << l2_error << std::endl;
        std::cout << "  L∞ norm:  " << linf_error << std::endl;
        
        // 2. Get error data for custom analysis
        std::cout << "\n2. Computing per-cell error data..." << std::endl;
        std::vector<Vector3D> gradDiffs = visualizer.computeGradientDifferences();
        std::vector<double> errorNorms = visualizer.computeErrorNormsPerCell();
        
        std::cout << "Computed " << gradDiffs.size() << " gradient differences" << std::endl;
        std::cout << "Computed " << errorNorms.size() << " error norms" << std::endl;
        
        // 3. Add errors to mesh for visualization
        std::cout << "\n3. Adding error data to mesh..." << std::endl;
        visualizer.add_errors_in_mesh("GradError", "ErrorMagnitude");
        
        // Save enhanced mesh
        visualizer.saveVisualization("complete_error_analysis.vtp");
        
        // 4. Generate epsilon analysis plot
        std::cout << "\n4. Generating epsilon analysis plot..." << std::endl;
        PlotGenerator generator("plots");
        
        std::string plotPath = generator.generateEpsilonErrorPlot(
            mesh, f_2, grad_f_2, kernel_4,
            0.01, 1.5, 25, NormType::L2,
            "f true grad", "Grad"
        );
        
        std::cout << "Generated epsilon plot: " << plotPath << std::endl;
        
        std::cout << "\n=== Analysis Complete ===" << std::endl;
        std::cout << "Files generated:" << std::endl;
        std::cout << "  • complete_error_analysis.vtp - Mesh with error data" << std::endl;
        std::cout << "  • " << plotPath << " - Epsilon analysis plot" << std::endl;
        
    } catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return EXIT_FAILURE;
    }
    
    return EXIT_SUCCESS;
}