#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <cmath>

#include "vector_3d.h"

// f_1(x) = x[0]
double f_1(const Vector3D &vec) { return vec.x; }

Vector3D grad_f_1(const Vector3D &vec) { return {1.0, 0, 0}; }

// f_2(x) = cos(\theta)
double f_2(const Vector3D &vec) {
    return -vec.x / std::sqrt(vec.x * vec.x + vec.z * vec.z);
}

Vector3D grad_f_2(const Vector3D &vec) {
    double denominator = std::pow(vec.x * vec.x + vec.z * vec.z, 1.5);
    return {-vec.z * vec.z / denominator, 0, vec.x * vec.z / denominator};
}

#endif