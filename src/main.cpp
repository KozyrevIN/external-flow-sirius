#include <cstdlib>
#include <iostream>
#include <string>


#include <vtkCellData.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkXMLPolyDataWriter.h>


#include "../include/functions.h"
#include "../include/geometry.h"
#include "../include/kernels.h"
#include "../include/utils.h"
#include "../include/visualizer.h"


const std::string mesh_in_path = "meshes/sphere_642.vtk";
const std::string mesh_out_path = "meshes/sphere_642_w_fields.vtp";

int main(int argc, char *argv[]) {
    vtkSmartPointer<vtkPolyData> mesh = load_mesh(mesh_in_path, false);

    // Create necessary attributes
    attach_center_to_cells(mesh);
    attach_area(mesh);

    // Compute true gradient and attach function values
    auto true_grad_array = compute_f_true_grad(mesh, grad_f_1);
    mesh->GetCellData()->AddArray(true_grad_array);
    attach_f(mesh, f_1);

    // Compute approximated gradient using grad_calculator
    grad_calculator grad_calc(f_1, 0.1, kernel_2);
    grad_calc.attach_grad(mesh);

    // Create and use the visualizer
    try {
        GradientVisualizer visualizer(mesh);

        // Compute gradient differences
        visualizer.computeGradientDifference("GradientError");

        // Compute per-cell error norms
        visualizer.computeErrorNormPerCell(NormType::L1, "ErrorNorm_L1");
        visualizer.computeErrorNormPerCell(NormType::L2, "ErrorNorm_L2");
        visualizer.computeErrorNormPerCell(NormType::LINF, "ErrorNorm_Linf");

        // Print error statistics
        visualizer.printErrorStatistics();

        // Save visualization with all error analysis
        visualizer.saveVisualization(mesh_out_path);

    } catch (const std::exception &e) {
        std::cerr << "Visualizer error: " << e.what() << std::endl;
    }
    return EXIT_SUCCESS;
}
