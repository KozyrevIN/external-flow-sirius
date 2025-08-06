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