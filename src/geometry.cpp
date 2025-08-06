#include "../include/geometry.h"
#include "../include/utils.h"
#include <vtkTriangle.h>
#include <vtkQuad.h>
#include <vtkPolygon.h>
#include <vtkMath.h>
#include <vtkPoints.h>
#include <iostream>

void get_cell_center(vtkCell *cell, std::vector<double> *center) {
    cell->GetParametricCenter(center->data());
}

void add_projection(vtkCell *cell, std::vector<double> *point, vtkPolyData* polyData, vtkIdType cellId) {
    vtkDataArray* projArray = get_cell_attribute(polyData, "projection");
    if (!projArray) {
        add_cell_attribute("projection", polyData, 3, VTK_DOUBLE);
        projArray = get_cell_attribute(polyData, "projection");
    }
    
    // Вычисляем нормаль к плоскости ячейки
    double normal[3] = {0.0, 0.0, 0.0};
    
    if (cell->GetCellType() == VTK_TRIANGLE) {
        // Для треугольника: нормаль через векторное произведение двух сторон
        vtkPoints* points = cell->GetPoints();
        if (points->GetNumberOfPoints() >= 3) {
            double p0[3], p1[3], p2[3];
            points->GetPoint(0, p0);
            points->GetPoint(1, p1);
            points->GetPoint(2, p2);
            
            // Векторы сторон треугольника
            double v1[3] = {p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]};
            double v2[3] = {p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]};
            
            // Векторное произведение дает нормаль
            vtkMath::Cross(v1, v2, normal);
            vtkMath::Normalize(normal);
        }
    }
    
    // Вычисляем проекцию вектора на плоскость ячейки
    // Проекция = вектор - (вектор · нормаль) * нормаль
    double pointVec[3] = {(*point)[0], (*point)[1], (*point)[2]};
    double dotProduct = vtkMath::Dot(pointVec, normal);
    
    double projection[3];
    projection[0] = pointVec[0] - dotProduct * normal[0];
    projection[1] = pointVec[1] - dotProduct * normal[1];
    projection[2] = pointVec[2] - dotProduct * normal[2];
    
    projArray->SetTuple3(cellId, projection[0], projection[1], projection[2]);
}