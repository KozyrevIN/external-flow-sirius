#include <cstdlib>
#include <iostream>

#include <vtkCellData.h>

#include "../include/functions.h"
#include "../include/geometry.h"
#include "../include/kernels.h"
#include "../include/utils.h"
#include "../include/visualizer.h"

const std::string mesh_in_path = "meshes/sphere_642.vtk";

int main(int argc, char *argv[]) {
    vtkSmartPointer<vtkPolyData> mesh = load_mesh(mesh_in_path, false);

    // Create necessary attributes
    attach_center_to_cells(mesh);
    attach_area(mesh);

    // Compute true gradient and attach function values
    auto true_grad_array = compute_f_true_grad(mesh, grad_f_1);
    mesh->GetCellData()->AddArray(true_grad_array);
    attach_f(mesh, f_1);

    // Compute gradient with grad_calculator to have "Grad" array
    grad_calculator grad_calc(f_1, 0.25, kernel_2);
    grad_calc.attach_grad(mesh);

    try {
        // Create visualizer
        GradientVisualizer visualizer(mesh);

        // Evaluate epsilon error for a range of epsilon values
        std::cout << "Starting epsilon error evaluation..." << std::endl;
        
        visualizer.evaluateEpsilonError(
            f_1,           // function
            kernel_2,      // kernel function
            0.01,          // epsilon_min
            1.0,           // epsilon_max
            20,            // num_points
            "epsilon_error_analysis.csv",  // csv filename
            NormType::L2   // norm type
        );

        std::cout << "Epsilon error evaluation completed!" << std::endl;

    } catch (const std::exception &e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}