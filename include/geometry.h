#include <vtkCell.h>
#include <vtkPolyData.h>
#include <functional>
#include "vector_3d.h"

void get_cell_center(vtkCell *cell, Vector3D *center);
void attach_center_to_cells(vtkPolyData* polyData);
Vector3D calc_projection(vtkCell *cell, Vector3D *vector);
Vector3D calc_normal(vtkCell *cell);
double calc_area(vtkCell *cell);
void attach_area(vtkPolyData* polyData);
Vector3D calc_grad(vtkPolyData* polyData, vtkIdType cellId, std::function<double(Vector3D)> f);