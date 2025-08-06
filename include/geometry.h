#include <vtkCell.h>
#include <vtkPolyData.h>
#include "vector_3d.h"

void get_cell_center(vtkCell *cell, Vector3D *center);
void add_projection(vtkCell *cell, Vector3D *point, vtkPolyData* polyData, vtkIdType cellId);
void add_center_to_cells(vtkPolyData* polyData);