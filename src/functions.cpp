#include <vtkCell.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>

#include "../include/functions.h"
#include "../include/geometry.h"
#include <cmath>
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

double f_3(const Vector3D &vec) {
    return vec.x * vec.x + vec.y * vec.y + vec.z * vec.z;
}

Vector3D grad_f_3(const Vector3D &vec) {
    return {2 * vec.x, 2 * vec.y, 2 * vec.z};
}

double f_4(const Vector3D &vec) { return exp(vec.y * vec.y - vec.z * vec.z); }

Vector3D grad_f_4(const Vector3D &vec) {
    return {0, 2 * vec.y * f_4(vec), -2 * vec.z * f_4(vec)};
}

double f_5(const Vector3D &vec) {
    double x_1 = 3.0;
    double k_1 = 1.5;
    double k_2 = 2.0;
    double k_3 = -1;
    double dist = std::sqrt((vec.x - x_1) * (vec.x - x_1) + vec.y * vec.y +
                            vec.z * vec.z);
    return cos((k_1 * (vec.x - x_1) + k_2 * vec.y + k_3 * vec.z)) / dist;
}

Vector3D grad_f_5(const Vector3D &vec) {
    double x_1 = 3.0;
    double k_1 = 1.5;
    double k_2 = 2.0;
    double k_3 = -1;

    const double r_x = vec.x - x_1;
    const double r_y = vec.y;
    const double r_z = vec.z;

    // 2. Calculate the squared magnitude of r, which is the common base of the denominator.
    // r_squared = (x - x_1)^2 + y^2 + z^2
    const double r_squared = r_x * r_x + r_y * r_y + r_z * r_z;

    // --- Safety Check: Avoid division by zero ---
    // If the distance is zero, the gradient is undefined. Return a zero vector.
    if (r_squared == 0.0) {
        return {0.0, 0.0, 0.0};
    }

    // 3. Calculate the argument for sin and cos, which is the dot product of k and r.
    // arg = k_1*(x-x_1) + k_2*y + k_3*z
    const double arg = k_1 * r_x + k_2 * r_y + k_3 * r_z;

    // 4. Calculate the full denominator, which is r^3.
    // Using pow is clear, but sqrt(r_squared) * r_squared could also be used.
    const double denominator = pow(r_squared, 1.5);

    // 5. Calculate sin(arg) and cos(arg) once.
    const double sin_arg = sin(arg);
    const double cos_arg = cos(arg);

    // 6. Calculate each component of the gradient using the pre-calculated terms.
    Vector3D gradient;

    gradient.x = -(k_1 * r_squared * sin_arg + r_x * cos_arg) / denominator;
    gradient.y = -(k_2 * r_squared * sin_arg + r_y * cos_arg) / denominator;
    gradient.z = -(k_3 * r_squared * sin_arg + r_z * cos_arg) / denominator;

    return gradient;
}

double f_5(const Vector3D &vec) {
    const Vector3D k = {1,1,1};
    const Vector3D q = {1,1,1};
    return cos(k.x * vec.x + k.y * vec.y + k.z * vec.z);
}
Vector3D grad_f_5(const Vector3D &vec) {
    return {0,0,0};
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

void attach_f_true_grad(vtkSmartPointer<vtkPolyData> mesh,
                        const std::function<Vector3D(Vector3D)> &f_grad) {
    vtkSmartPointer<vtkDoubleArray> true_grad_array =
        compute_f_true_grad(mesh, f_grad);
    mesh->GetCellData()->AddArray(true_grad_array);
}
