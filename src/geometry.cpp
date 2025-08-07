#include "../include/geometry.h"
#include "../include/utils.h"
#include <vtkTriangle.h>
#include <vtkQuad.h>
#include <vtkPolygon.h>
#include <vtkMath.h>
#include <vtkPoints.h>
#include <vtkFloatArray.h>
#include <vtkSmartPointer.h>
#include <vtkCellData.h>
#include <iostream>
#include <cmath>
#include <functional>
#include "../include/vector_3d.h"
#include "../include/utils.h"

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
void attach_center_to_cells(vtkPolyData* polyData) {
    vtkIdType numCells = polyData->GetNumberOfCells();
    
    // Создаем массив для хранения центров ячеек
    vtkSmartPointer<vtkDoubleArray> centerArray = vtkSmartPointer<vtkDoubleArray>::New();
    centerArray->SetNumberOfComponents(3);
    centerArray->SetNumberOfTuples(numCells);
    centerArray->SetName("CellCenters");
    
    // Вычисляем центр для каждой ячейки
    for (vtkIdType i = 0; i < numCells; i++) {
        vtkCell* cell = polyData->GetCell(i);
        Vector3D center;
        get_cell_center(cell, &center);
        
        // Сохраняем координаты центра в массив
        double centerCoords[3] = {
            static_cast<double>(center.x),
            static_cast<double>(center.y), 
            static_cast<double>(center.z)
        };
        centerArray->SetTuple(i, centerCoords);
    }
    
    // Добавляем массив как атрибут к данным ячеек
    polyData->GetCellData()->AddArray(centerArray);
    
    // Отладочная информация
    std::cout << "Добавлен массив центров ячеек: " << numCells << " центров" << std::endl;
    std::cout << "Количество массивов в CellData: " << polyData->GetCellData()->GetNumberOfArrays() << std::endl;
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
Vector3D calc_projection(vtkCell *cell, Vector3D *vector) {
    // Получаем точки треугольника
    double p0[3], p1[3], p2[3];
    cell->GetPoints()->GetPoint(0, p0);
    cell->GetPoints()->GetPoint(1, p1);
    cell->GetPoints()->GetPoint(2, p2);
    
    // Вычисляем нормаль треугольника через метод vtkTriangle
    double normal[3];
    vtkTriangle::ComputeNormal(p0, p1, p2, normal);
    
    // Вычисляем скалярное произведение вектора на нормаль
    double dot_product = vector->x * normal[0] + vector->y * normal[1] + vector->z * normal[2];
    
    // Проекция вектора на плоскость = вектор - (нормальная составляющая)
    Vector3D projection;
    projection.x = vector->x - dot_product * normal[0];
    projection.y = vector->y - dot_product * normal[1];
    projection.z = vector->z - dot_product * normal[2];
    
    return projection;
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

void attach_area(vtkPolyData* polyData) {
   
    vtkIdType numCells = polyData->GetNumberOfCells();
    
    // Создаем массив для хранения площадей ячеек
    vtkSmartPointer<vtkDoubleArray> areaArray = vtkSmartPointer<vtkDoubleArray>::New();
    areaArray->SetNumberOfComponents(1);
    areaArray->SetNumberOfTuples(numCells);
    areaArray->SetName("Area");
    
    // Вычисляем площадь для каждой ячейки
    for (vtkIdType i = 0; i < numCells; i++) {
        vtkCell* cell = polyData->GetCell(i);
        double area = calc_area(cell);
        areaArray->SetTuple1(i, area);
    }
    
    // Добавляем массив как атрибут к данным ячеек
    polyData->GetCellData()->AddArray(areaArray);
    
    // Отладочная информация
}

grad_calculator::grad_calculator(std::function<double(Vector3D)> f, double epsilon, std::function<double(double)> kernel) : f(f), epsilon(epsilon), kernel(kernel) {}

Vector3D grad_calculator::calc_grad(vtkSmartPointer<vtkPolyData> polyData, vtkIdType cellId) {
        vtkCell* cell = polyData->GetCell(cellId);
        Vector3D grad(0.0, 0.0, 0.0);
        double area = getAttributeArea(cellId, polyData);
        Vector3D center = getCenter(cellId, polyData);
        double f_value = f(center);
        for (vtkIdType i = 0; i < cell->GetNumberOfPoints(); i++) {
            Vector3D center_i = getCenter(i, polyData);
            double area_i = getAttributeArea(i, polyData);
            double f_value_i = f(center_i);
            double distance = sqrt(pow(center_i.x - center.x, 2) + pow(center_i.y - center.y, 2) + pow(center_i.z - center.z, 2));
            double coeff = (f_value_i - f_value) * area_i * kernel(distance);
            Vector3D diff(center_i.x - center.x, center_i.y - center.y, center_i.z - center.z);
            Vector3D projection = calc_projection(cell, &diff);
            grad.x += coeff * projection.x;
            grad.y += coeff * projection.y;
            grad.z += coeff * projection.z;
        }
        grad.x /= pow(epsilon,4);
        grad.y /= pow(epsilon,4);
        grad.z /= pow(epsilon,4);
        
        return grad;
}

void grad_calculator::attach_grad(vtkSmartPointer<vtkPolyData> polyData) {
        vtkIdType numCells = polyData->GetNumberOfCells();
        vtkSmartPointer<vtkDoubleArray> gradArray = vtkSmartPointer<vtkDoubleArray>::New();
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