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
void add_center_to_cells(vtkPolyData* polyData) {
    vtkIdType numCells = polyData->GetNumberOfCells();
    
    // Создаем массив для хранения центров ячеек
    vtkSmartPointer<vtkFloatArray> centerArray = vtkSmartPointer<vtkFloatArray>::New();
    centerArray->SetNumberOfComponents(3);
    centerArray->SetNumberOfTuples(numCells);
    centerArray->SetName("CellCenters");
    
    // Вычисляем центр для каждой ячейки
    for (vtkIdType i = 0; i < numCells; i++) {
        vtkCell* cell = polyData->GetCell(i);
        Vector3D center;
        get_cell_center(cell, &center);
        
        // Сохраняем координаты центра в массив
        float centerCoords[3] = {
            static_cast<float>(center.x),
            static_cast<float>(center.y), 
            static_cast<float>(center.z)
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
    vtkSmartPointer<vtkFloatArray> areaArray = vtkSmartPointer<vtkFloatArray>::New();
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
    std::cout << "Добавлен массив площадей ячеек: " << numCells << " площадей" << std::endl;
}

Vector3D calc_grad(vtkPolyData* polyData, vtkIdType cellId, std::function<double(Vector3D)> f) {
    vtkCell* cell = polyData->GetCell(cellId);
    Vector3D grad(0.0, 0.0, 0.0);
    
    for (vtkIdType i = 0; i < cell->GetNumberOfPoints(); i++) {
        // Получаем координаты точки
        double coords[3];
        cell->GetPoints()->GetPoint(i, coords);
        Vector3D point(coords[0], coords[1], coords[2]);
        
        // Вычисляем градиент
        Vector3D normal = calc_normal(cell);
        double f_value = f(point);
        
        grad.x += f_value * normal.x;
        grad.y += f_value * normal.y;
        grad.z += f_value * normal.z;
    }
    
    // Нормируем на площадь
    double area = calc_area(cell);
    if (area > 0) {
        grad.x /= area;
        grad.y /= area;
        grad.z /= area;
    }
    
    return grad;
}