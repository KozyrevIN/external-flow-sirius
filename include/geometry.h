#include <vtkCell.h>
#include <vtkPolyData.h>
#include <vector>

void get_cell_center(vtkCell *cell, std::vector<double> *center);

// Вычисляет проекцию вектора point на плоскость ячейки и добавляет как атрибут "projection"
void add_projection(vtkCell *cell, std::vector<double> *point, vtkPolyData* polyData, vtkIdType cellId);