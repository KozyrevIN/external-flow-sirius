#include "../include/visualizer.h"
#include "../include/geometry.h"
#include "../include/utils.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>
#include <vtkXMLPolyDataWriter.h>


GradientVisualizer::GradientVisualizer(vtkSmartPointer<vtkPolyData> mesh,
                                       const std::string &trueGradArrayName,
                                       const std::string &computedGradArrayName)
    : mesh(mesh), true_grad_array_name(trueGradArrayName),
      computed_grad_array_name(computedGradArrayName) {

    if (!mesh) {
        throw std::runtime_error("Mesh cannot be null");
    }

    // Verify that the required arrays exist
    if (!mesh->GetCellData()->GetArray(trueGradArrayName.c_str())) {
        throw std::runtime_error("True gradient array '" + trueGradArrayName +
                                 "' not found in mesh");
    }

    if (!mesh->GetCellData()->GetArray(computedGradArrayName.c_str())) {
        throw std::runtime_error("Computed gradient array '" +
                                 computedGradArrayName + "' not found in mesh");
    }
}

Vector3D GradientVisualizer::getVectorFromArray(vtkDoubleArray *array,
                                                vtkIdType cellId) {
    double vec[3];
    array->GetTuple(cellId, vec);
    return Vector3D(vec[0], vec[1], vec[2]);
}

void GradientVisualizer::setVectorToArray(vtkDoubleArray *array,
                                          vtkIdType cellId,
                                          const Vector3D &vec) {
    double vecArray[3] = {vec.x, vec.y, vec.z};
    array->SetTuple(cellId, vecArray);
}

double GradientVisualizer::computeVectorNorm(const Vector3D &vec,
                                             NormType normType) {
    // Individual vector norm is always Euclidean (L2)
    return std::sqrt(vec.x * vec.x + vec.y * vec.y + vec.z * vec.z);
}

std::vector<Vector3D> GradientVisualizer::computeGradientDifferences() {
    vtkDoubleArray *trueGradArray = vtkDoubleArray::SafeDownCast(
        mesh->GetCellData()->GetArray(true_grad_array_name.c_str()));
    vtkDoubleArray *computedGradArray = vtkDoubleArray::SafeDownCast(
        mesh->GetCellData()->GetArray(computed_grad_array_name.c_str()));

    if (!trueGradArray || !computedGradArray) {
        throw std::runtime_error("Failed to retrieve gradient arrays");
    }

    vtkIdType numCells = mesh->GetNumberOfCells();
    std::vector<Vector3D> differences;
    differences.reserve(numCells);

    // Compute differences for each cell
    for (vtkIdType i = 0; i < numCells; i++) {
        Vector3D trueGrad = getVectorFromArray(trueGradArray, i);
        Vector3D computedGrad = getVectorFromArray(computedGradArray, i);

        differences.emplace_back(
            trueGrad.x - computedGrad.x,
            trueGrad.y - computedGrad.y,
            trueGrad.z - computedGrad.z
        );
    }

    return differences;
}

double GradientVisualizer::computeErrorNorm(NormType normType) {
    // Check if area data exists
    vtkDataArray *areaArray = mesh->GetCellData()->GetArray("Area");
    if (!areaArray) {
        throw std::runtime_error("Failed to retrieve area array. Ensure mesh has area data attached.");
    }

    // Use computeErrorNormsPerCell for consistency and reuse
    std::vector<double> cellErrors = computeErrorNormsPerCell();
    vtkIdType numCells = mesh->GetNumberOfCells();
    
    double result = 0.0;

    switch (normType) {
    case NormType::L1: {
        // L1 norm: ∫|e(x)| dA
        for (vtkIdType i = 0; i < numCells; i++) {
            double area = areaArray->GetTuple1(i);
            result += cellErrors[i] * area;
        }
        return result;
    }
    case NormType::L2: {
        // L2 norm: sqrt(∫|e(x)|² dA)
        for (vtkIdType i = 0; i < numCells; i++) {
            double area = areaArray->GetTuple1(i);
            result += cellErrors[i] * cellErrors[i] * area;
        }
        return std::sqrt(result);
    }
    case NormType::LINF: {
        // L∞ norm is not area-weighted (maximum over all cells)
        for (vtkIdType i = 0; i < numCells; i++) {
            result = std::max(result, cellErrors[i]);
        }
        return result;
    }
    default:
        throw std::runtime_error("Unknown norm type");
    }
}

std::vector<double> GradientVisualizer::computeErrorNormsPerCell() {
    vtkDoubleArray *trueGradArray = vtkDoubleArray::SafeDownCast(
        mesh->GetCellData()->GetArray(true_grad_array_name.c_str()));
    vtkDoubleArray *computedGradArray = vtkDoubleArray::SafeDownCast(
        mesh->GetCellData()->GetArray(computed_grad_array_name.c_str()));

    if (!trueGradArray || !computedGradArray) {
        throw std::runtime_error("Failed to retrieve gradient arrays");
    }

    vtkIdType numCells = mesh->GetNumberOfCells();
    std::vector<double> errorNorms;
    errorNorms.reserve(numCells);

    for (vtkIdType i = 0; i < numCells; i++) {
        Vector3D trueGrad = getVectorFromArray(trueGradArray, i);
        Vector3D computedGrad = getVectorFromArray(computedGradArray, i);

        double dx = trueGrad.x - computedGrad.x;
        double dy = trueGrad.y - computedGrad.y;
        double dz = trueGrad.z - computedGrad.z;

        // Individual vector norm is always Euclidean (L2)
        errorNorms.push_back(std::sqrt(dx*dx + dy*dy + dz*dz));
    }

    return errorNorms;
}

void GradientVisualizer::saveVisualization(const std::string &filename) {
    vtkSmartPointer<vtkXMLPolyDataWriter> writer =
        vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer->SetFileName(filename.c_str());
    writer->SetInputData(mesh);
    writer->Write();

    std::cout << "Saved visualization to: " << filename << std::endl;
}

void GradientVisualizer::printErrorStatistics() {
    try {
        double l1Error = computeErrorNorm(NormType::L1);
        double l2Error = computeErrorNorm(NormType::L2);
        double linfError = computeErrorNorm(NormType::LINF);

        std::cout << "\n=== Gradient Error Statistics ===" << std::endl;
        std::cout << "L1 norm:  " << l1Error << std::endl;
        std::cout << "L2 norm:  " << l2Error << std::endl;
        std::cout << "L∞ norm:  " << linfError << std::endl;

        vtkIdType numCells = mesh->GetNumberOfCells();
        std::cout << "Average L1 error per cell: " << l1Error / numCells
                  << std::endl;
        std::cout << "Average L2 error per cell: " << l2Error / numCells
                  << std::endl;
        std::cout << "==================================" << std::endl;

    } catch (const std::exception &e) {
        std::cerr << "Error computing statistics: " << e.what() << std::endl;
    }
}

void GradientVisualizer::setTrueGradientArrayName(const std::string &name) {
    true_grad_array_name = name;
}

void GradientVisualizer::setComputedGradientArrayName(const std::string &name) {
    computed_grad_array_name = name;
}

void GradientVisualizer::evaluateEpsilonError(
    std::function<double(Vector3D)> function,
    std::function<double(double)> kernel,
    double epsilon_min,
    double epsilon_max,
    int num_points,
    const std::string& csv_filename,
    NormType normType) {
    
    std::ofstream csv_file(csv_filename);
    if (!csv_file.is_open()) {
        throw std::runtime_error("Could not open CSV file: " + csv_filename);
    }
    
    // Write CSV header
    std::string norm_name;
    switch (normType) {
        case NormType::L1: norm_name = "L1"; break;
        case NormType::L2: norm_name = "L2"; break;
        case NormType::LINF: norm_name = "Linf"; break;
    }
    csv_file << "epsilon," << norm_name << "_error\n";
    
    std::cout << "Evaluating gradient error for epsilon range [" 
              << epsilon_min << ", " << epsilon_max << "] with " 
              << num_points << " points..." << std::endl;
    
    // Store original computed gradient array name to restore later
    std::string original_computed_name = computed_grad_array_name;
    
    for (int i = 0; i < num_points; i++) {
        // Calculate epsilon value (logarithmic spacing for better coverage)
        double epsilon = epsilon_min * std::pow(epsilon_max / epsilon_min, 
                                               static_cast<double>(i) / (num_points - 1));
        
        // Create grad_calculator with current epsilon
        grad_calculator grad_calc(function, epsilon, kernel);
        
        // Create a temporary mesh copy to avoid modifying the original
        vtkSmartPointer<vtkPolyData> temp_mesh = vtkSmartPointer<vtkPolyData>::New();
        temp_mesh->DeepCopy(mesh);
        
        // Compute gradient with current epsilon
        grad_calc.attach_grad(temp_mesh);
        
        // Temporarily update computed gradient array name for error calculation
        std::string temp_computed_name = "Grad";
        computed_grad_array_name = temp_computed_name;
        
        // Create temporary visualizer for this epsilon
        GradientVisualizer temp_visualizer(temp_mesh, true_grad_array_name, temp_computed_name);
        
        // Compute error norm
        double error = temp_visualizer.computeErrorNorm(normType);
        
        // Write to CSV
        csv_file << std::fixed << std::setprecision(6) << epsilon << "," << error << "\n";
    }
    
    // Restore original computed gradient array name
    computed_grad_array_name = original_computed_name;
    
    csv_file.close();
    std::cout << "Epsilon error analysis saved to: " << csv_filename << std::endl;
}

void GradientVisualizer::add_errors_in_mesh(const std::string &differenceArrayName,
                                           const std::string &errorNormArrayName) {
    // Ensure meshes directory exists
    std::string mkdir_cmd = "mkdir -p meshes";
    std::system(mkdir_cmd.c_str());
    
    // Compute gradient differences and error norms
    std::vector<Vector3D> differences = computeGradientDifferences();
    std::vector<double> errorNorms = computeErrorNormsPerCell();
    
    vtkIdType numCells = mesh->GetNumberOfCells();
    
    // Create difference array (vector)
    vtkSmartPointer<vtkDoubleArray> diffArray = vtkSmartPointer<vtkDoubleArray>::New();
    diffArray->SetNumberOfComponents(3);
    diffArray->SetNumberOfTuples(numCells);
    diffArray->SetName(differenceArrayName.c_str());
    
    // Create error norm array (scalar)
    vtkSmartPointer<vtkDoubleArray> errorArray = vtkSmartPointer<vtkDoubleArray>::New();
    errorArray->SetNumberOfComponents(1);
    errorArray->SetNumberOfTuples(numCells);
    errorArray->SetName(errorNormArrayName.c_str());
    
    // Fill arrays
    for (vtkIdType i = 0; i < numCells; i++) {
        // Set difference vector
        double diffVec[3] = {differences[i].x, differences[i].y, differences[i].z};
        diffArray->SetTuple(i, diffVec);
        
        // Set error norm scalar
        errorArray->SetTuple1(i, errorNorms[i]);
    }
    
    // Add arrays to mesh
    mesh->GetCellData()->AddArray(diffArray);
    mesh->GetCellData()->AddArray(errorArray);
    
    std::cout << "Added gradient differences and error norms to mesh for " 
              << numCells << " cells" << std::endl;
}