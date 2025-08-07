#include "vector_3d.h"
#include <functional>
#include <iostream>
#include <vtkCell.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>

void get_cell_center(vtkCell *cell, Vector3D *center);
void attach_center_to_cells(vtkPolyData *polyData);
Vector3D calc_projection(vtkCell *cell, Vector3D *vector);
Vector3D calc_normal(vtkCell *cell);
double calc_area(vtkCell *cell);
void attach_area(vtkPolyData* polyData);

struct grad_calculator{
    std::function<double(Vector3D)> f;
    double epsilon;
    std::function<double(double)> kernel;
    grad_calculator(std::function<double(Vector3D)> f, double epsilon, std::function<double(double)> kernel);
    Vector3D calc_grad(vtkSmartPointer<vtkPolyData> polyData, vtkIdType cellId);
    void attach_grad(vtkSmartPointer<vtkPolyData> polyData);
};