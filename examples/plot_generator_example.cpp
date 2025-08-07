#include <cstdlib>
#include <iostream>

#include "../include/functions.h"
#include "../include/geometry.h"
#include "../include/kernels.h"
#include "../include/plot_generator.h"
#include "../include/utils.h"

const std::string mesh_in_path = "meshes/sphere_642.vtk";


int main(int argc, char *argv[]) {
    try {
        // Load and prepare mesh
        vtkSmartPointer<vtkPolyData> mesh = load_mesh(mesh_in_path, false);
        attach_center_to_cells(mesh);
        attach_area(mesh);
        
        // Attach true gradient and function values
        attach_f_true_grad(mesh, grad_f_2);
        attach_f(mesh, f_2);
        
        // Calculate computed gradient with some epsilon
        grad_calculator grad_calc(f_2, 0.25, kernel_2);
        grad_calc.attach_grad(mesh);
        
        std::cout << "=== Testing PlotGenerator ===" << std::endl;
        
        // Method 1: Using PlotGenerator class
        PlotGenerator generator("plots", "plot_epsilon_error.py");
        
        std::string plot1 = generator.generateEpsilonErrorPlot(
            mesh, f_2, grad_f_2, kernel_2,
            0.01, 1.0, 15, NormType::L2
        );
        
        std::cout << "Generated plot 1: " << plot1 << std::endl;
        
        // Method 2: Using static quickPlot method
        std::string plot2 = PlotGenerator::quickPlot(
            mesh, f_2, grad_f_2, kernel_4,
            0.05, 0.8
        );
        
        std::cout << "Generated plot 2: " << plot2 << std::endl;
        
        // Method 3: Different norm type
        std::string plot3 = generator.generateEpsilonErrorPlot(
            mesh, f_2, grad_f_2, kernel_6,
            0.01, 2.0, 25, NormType::L1
        );
        
        std::cout << "Generated plot 3: " << plot3 << std::endl;
        

        std::cout << "\n=== All plots generated successfully! ===" << std::endl;
        std::cout << "Check the 'plots/' directory for PNG files." << std::endl;
        
    } catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return EXIT_FAILURE;
    }
    
    return EXIT_SUCCESS;
}