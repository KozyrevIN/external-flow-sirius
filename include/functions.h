#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <cmath>
#include <functional>

#include <vtkPolyData.h>
#include <vtkSmartPointer.h>

#include "vector_3d.h"

// f_1(x) = x[0]
double f_1(const Vector3D &vec);

Vector3D grad_f_1(const Vector3D &vec);

// f_2(x) = cos(\theta)
double f_2(const Vector3D &vec);

Vector3D grad_f_2(const Vector3D &vec);

// functions to attach f anf grad f to mesh

void attach_f(vtkSmartPointer<vtkPolyData> mesh,
              const std::function<double(Vector3D)> &f);

vtkSmartPointer<vtkDoubleArray> compute_f_true_grad(vtkSmartPointer<vtkPolyData> mesh,
                                                   const std::function<Vector3D(Vector3D)> &f_grad);

#endif