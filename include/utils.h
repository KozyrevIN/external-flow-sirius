#ifndef UTILS_H
#define UTILS_H

#include <string>
#include <vector>
#include <stdexcept>
#include <functional>

#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkDataArray.h>
#include <vtkFloatArray.h>
#include <vtkDoubleArray.h>
#include "vector_3d.h"

vtkSmartPointer<vtkPolyData> load_mesh(std::string filepath,
                                       bool verbose = false);
Vector3D getCenter(vtkIdType cellId, vtkPolyData* polyData);
double getAttributeArea(vtkIdType cellId, vtkPolyData* polyData);
vtkSmartPointer<vtkPolyData> load_and_init_mash(std::string mesh_in_path);
void write_mesh(vtkSmartPointer<vtkPolyData> mesh, std::string mesh_out_path);
void add_grads(vtkSmartPointer<vtkPolyData> mesh, std::function<double(Vector3D)> f, std::function<Vector3D(Vector3D)> grad_f, double epsilon, std::function<double(double)> kernel);

#endif