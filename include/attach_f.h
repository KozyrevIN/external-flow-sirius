#ifndef ATTACH_F_H
#define ATTACH_F_H

#include <functional>

#include <vtkPolyData.h>
#include <vtkSmartPointer.h>

#include "vector_3d.h"

void attach_f(vtkSmartPointer<vtkPolyData> mesh,
              const std::function<double(Vector3D)> &f);

#endif