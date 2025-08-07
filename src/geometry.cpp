#include "../include/geometry.h"
#include "../include/utils.h"
#include "../include/vector_3d.h"
#include <cmath>
#include <functional>
#include <iostream>
#include <vtkCellData.h>
#include <vtkFloatArray.h>
#include <vtkMath.h>
#include <vtkPoints.h>
#include <vtkPolygon.h>
#include <vtkQuad.h>
#include <vtkSmartPointer.h>
#include <vtkTriangle.h>

void get_cell_center(vtkCell *cell, Vector3D *center) {
    double coords[3];
    double pcoords[3] = {0.5, 0.5, 0.5}; // parametric coordinates for center
    double *weights = new double[cell->GetNumberOfPoints()];
    int subId = 0;

    cell->EvaluateLocation(subId, pcoords, coords, weights);

    center->x = coords[0];
    center->y = coords[1];
    center->z = coords[2];

    delete[] weights;
}
void attach_center_to_cells(vtkPolyData *polyData) {
    vtkIdType numCells = polyData->GetNumberOfCells();

    // Создаем массив для хранения центров ячеек
    vtkSmartPointer<vtkDoubleArray> centerArray =
        vtkSmartPointer<vtkDoubleArray>::New();
    centerArray->SetNumberOfComponents(3);
    centerArray->SetNumberOfTuples(numCells);
    centerArray->SetName("CellCenters");

    // Вычисляем центр для каждой ячейки
    for (vtkIdType i = 0; i < numCells; i++) {
        vtkCell *cell = polyData->GetCell(i);
        Vector3D center;
        get_cell_center(cell, &center);

        // Сохраняем координаты центра в массив
        double centerCoords[3] = {static_cast<double>(center.x),
                                  static_cast<double>(center.y),
                                  static_cast<double>(center.z)};
        centerArray->SetTuple(i, centerCoords);
    }

    // Добавляем массив как атрибут к данным ячеек
    polyData->GetCellData()->AddArray(centerArray);
}
Vector3D calc_normal(vtkCell *cell) {
    double p0[3], p1[3], p2[3];
    cell->GetPoints()->GetPoint(0, p0);
    cell->GetPoints()->GetPoint(1, p1);
    cell->GetPoints()->GetPoint(2, p2);

    // Вычисляем нормаль треугольника через метод vtkTriangle
    double normal[3];
    vtkTriangle::ComputeNormal(p0, p1, p2, normal);

    Vector3D normal_vector;
    normal_vector.x = normal[0];
    normal_vector.y = normal[1];
    normal_vector.z = normal[2];

    return normal_vector;
}
Vector3D calc_rectangle_projection(double* point1, double* point2, double* point3, Vector3D *vector) {
    // Получаем точки треугольника
    double p0[3] = {point1[0], point1[1], point1[2]};
    double p1[3] = {point2[0], point2[1], point2[2]};
    double p2[3] = {point3[0], point3[1], point3[2]};

    // Вычисляем нормаль треугольника через метод vtkTriangle
    double normal[3];
    vtkTriangle::ComputeNormal(p0, p1, p2, normal);

    // Вычисляем скалярное произведение вектора на нормаль
    double dot_product =
        vector->x * normal[0] + vector->y * normal[1] + vector->z * normal[2];

    // Проекция вектора на плоскость = вектор - (нормальная составляющая)
    Vector3D projection;
    projection.x = vector->x - dot_product * normal[0];
    projection.y = vector->y - dot_product * normal[1];
    projection.z = vector->z - dot_product * normal[2];

    return projection;
}
Vector3D calc_projection(vtkCell *cell, Vector3D *vector) {
    // Получаем точки треугольника
    double p0[3], p1[3], p2[3];
    cell->GetPoints()->GetPoint(0, p0);
    cell->GetPoints()->GetPoint(1, p1);
    cell->GetPoints()->GetPoint(2, p2);
    return calc_rectangle_projection(p0, p1, p2, vector);
}
double calc_area(vtkCell *cell) {
    // Получаем точки треугольника
    double p0[3], p1[3], p2[3];
    cell->GetPoints()->GetPoint(0, p0);
    cell->GetPoints()->GetPoint(1, p1);
    cell->GetPoints()->GetPoint(2, p2);

    // Используем встроенный метод VTK для вычисления площади треугольника
    return vtkTriangle::TriangleArea(p0, p1, p2);
}

void attach_area(vtkPolyData *polyData) {

    vtkIdType numCells = polyData->GetNumberOfCells();

    // Создаем массив для хранения площадей ячеек
    vtkSmartPointer<vtkDoubleArray> areaArray =
        vtkSmartPointer<vtkDoubleArray>::New();
    areaArray->SetNumberOfComponents(1);
    areaArray->SetNumberOfTuples(numCells);
    areaArray->SetName("Area");

    // Вычисляем площадь для каждой ячейки
    for (vtkIdType i = 0; i < numCells; i++) {
        vtkCell *cell = polyData->GetCell(i);
        double area = calc_area(cell);
        areaArray->SetTuple1(i, area);
    }

    // Добавляем массив как атрибут к данным ячеек
    polyData->GetCellData()->AddArray(areaArray);

    // Отладочная информация
}

grad_calculator::grad_calculator(std::function<double(Vector3D)> f,
                                 double epsilon,
                                 std::function<double(double)> kernel)
    : f(f), epsilon(epsilon), kernel(kernel) {}

Vector3D grad_calculator::calc_grad(vtkSmartPointer<vtkPolyData> polyData,
                                    vtkIdType cellId) const {
    vtkCell *cell = polyData->GetCell(cellId);
    Vector3D grad(0.0, 0.0, 0.0);
    Vector3D center = getCenter(cellId, polyData);
    double f_value = f(center);
    for (vtkIdType i = 0; i < polyData->GetNumberOfCells(); i++) {
        Vector3D center_i = getCenter(i, polyData);
        double area_i = getAttributeArea(i, polyData);
        double f_value_i = f(center_i);
        double distance =
            sqrt(pow(center_i.x - center.x, 2) + pow(center_i.y - center.y, 2) +
                 pow(center_i.z - center.z, 2));
        double coeff = (f_value_i - f_value) * area_i * kernel(distance/epsilon);
        Vector3D diff(center.x - center_i.x, center.y - center_i.y,
                      center.z - center_i.z);
        grad.x += coeff * diff.x;
        grad.y += coeff * diff.y;
        grad.z += coeff * diff.z;
    }
    grad = calc_projection(cell, &grad);
    grad.x /= pow(epsilon, 4);
    grad.y /= pow(epsilon, 4);
    grad.z /= pow(epsilon, 4);

    return grad;
}

void grad_calculator::attach_grad(vtkSmartPointer<vtkPolyData> polyData) {
    vtkIdType numCells = polyData->GetNumberOfCells();
    vtkSmartPointer<vtkDoubleArray> gradArray =
        vtkSmartPointer<vtkDoubleArray>::New();
    gradArray->SetNumberOfComponents(3);
    gradArray->SetNumberOfTuples(numCells);
    gradArray->SetName("Grad");
    for (vtkIdType i = 0; i < numCells; i++) {
        Vector3D grad = calc_grad(polyData, i);
        double gradValues[3] = {grad.x, grad.y, grad.z};
        gradArray->SetTuple(i, gradValues);
    }
    polyData->GetCellData()->SetVectors(gradArray);
}