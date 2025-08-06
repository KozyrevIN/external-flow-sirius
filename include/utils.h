#ifndef UTILS_H
#define UTILS_H

#include <string>
#include <vector>
#include <stdexcept>

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

#endif