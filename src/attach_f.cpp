#include <vtkCell.h>
#include <vtkCellData.h>
#include <vtkDoubleArray.h>

#include "../include/attach_f.h"

void attach_f(vtkSmartPointer<vtkPolyData> mesh,
              const std::function<double(Vector3D)> &f) {

    vtkSmartPointer<vtkDoubleArray> scalarArray =
        vtkSmartPointer<vtkDoubleArray>::New();
    scalarArray->SetName("f");
    scalarArray->SetNumberOfComponents(1);
    scalarArray->SetNumberOfTuples(mesh->GetNumberOfCells());

    for (vtkIdType cellId = 0; cellId < mesh->GetNumberOfCells(); ++cellId) {
        vtkCell* cell = mesh->GetCell(cellId);
        vtkIdList* pointIds = cell->GetPointIds();
        int numPointsInCell = pointIds->GetNumberOfIds();
        
        // Calculate the cell's center (barycenter)
        double center[3] = {0.0, 0.0, 0.0};
        for (vtkIdType i = 0; i < numPointsInCell; ++i) {
            double pointCoords[3];
            mesh->GetPoint(pointIds->GetId(i), pointCoords);
            center[0] += pointCoords[0];
            center[1] += pointCoords[1];
            center[2] += pointCoords[2];
        }
        center[0] /= numPointsInCell;
        center[1] /= numPointsInCell;
        center[2] /= numPointsInCell;

        // Evaluate our function f(x, y, z) at the center
        double value = f({center[0], center[1], center[2]});

        // Store this value in our array at the cell's index
        scalarArray->SetTuple1(cellId, value);
    }

    mesh->GetCellData()->AddArray(scalarArray);
}