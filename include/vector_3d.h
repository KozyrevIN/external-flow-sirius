#ifndef VECTOR_3D_H
#define VECTOR_3D_H

struct Vector3D {
    double x;
    double y;
    double z;
    Vector3D(double x, double y, double z) : x(x), y(y), z(z) {}
    Vector3D() : x(0), y(0), z(0) {}
    Vector3D(const Vector3D &other) : x(other.x), y(other.y), z(other.z) {}
};

#endif