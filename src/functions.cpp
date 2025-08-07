#include <vtkCell.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>

#include "../include/functions.h"
#include "../include/geometry.h"

// f_1(x) = x[0]
double f_1(const Vector3D &vec) { return vec.x; }

Vector3D grad_f_1(const Vector3D &vec) { return {1.0, 0, 0}; }

// f_2(x) = cos(\theta)
double f_2(const Vector3D &vec) {
    return -vec.x / std::sqrt(vec.x * vec.x + vec.y * vec.y + vec.z * vec.z);
}

Vector3D grad_f_2(const Vector3D &vec) {
    double denominator =
        std::pow(vec.x * vec.x + vec.y * vec.y + vec.z * vec.z, 1.5);
    return {-(vec.z * vec.z + vec.y * vec.y) / denominator,
            vec.x * vec.y / denominator, vec.x * vec.z / denominator};
}

void attach_f(vtkSmartPointer<vtkPolyData> mesh,
              const std::function<double(Vector3D)> &f) {

    vtkSmartPointer<vtkDoubleArray> f_vals =
        vtkSmartPointer<vtkDoubleArray>::New();
    f_vals->SetName("f");
    f_vals->SetNumberOfComponents(1);
    f_vals->SetNumberOfTuples(mesh->GetNumberOfCells());

    vtkSmartPointer<vtkDataArray> centers =
        mesh->GetCellData()->GetArray("CellCenters");

    for (vtkIdType cellId = 0; cellId < mesh->GetNumberOfCells(); ++cellId) {
        // Get the cell's center
        double center[3];
        centers->GetTuple(cellId, center);

        // Evaluate our function f(x, y, z) at the center
        double value = f({center[0], center[1], center[2]});

        // Store this value in our array at the cell's index
        f_vals->SetTuple1(cellId, value);
    }

    mesh->GetCellData()->AddArray(f_vals);
}

vtkSmartPointer<vtkDoubleArray>
compute_f_true_grad(vtkSmartPointer<vtkPolyData> mesh,
                    const std::function<Vector3D(Vector3D)> &f_grad) {

    vtkSmartPointer<vtkDoubleArray> vectorArray =
        vtkSmartPointer<vtkDoubleArray>::New();
    vectorArray->SetName("f true grad");
    vectorArray->SetNumberOfComponents(3);
    vectorArray->SetNumberOfTuples(mesh->GetNumberOfCells());
    vtkSmartPointer<vtkDataArray> centers =
        mesh->GetCellData()->GetArray("CellCenters");

    for (vtkIdType cellId = 0; cellId < mesh->GetNumberOfCells(); ++cellId) {
        // Get the cell's center from pre-defined CellCenters array
        double center[3];
        centers->GetTuple(cellId, center);

        // Evaluate our function f_grad(x, y, z) at the center
        Vector3D grad = f_grad({center[0], center[1], center[2]});
        Vector3D grad_projected = calc_projection(mesh->GetCell(cellId), &grad);

        // Store this value in our array at the cell's index
        double gradArray[3] = {grad_projected.x, grad_projected.y,
                               grad_projected.z};
        vectorArray->SetTuple(cellId, gradArray);
    }

    return vectorArray;
}
