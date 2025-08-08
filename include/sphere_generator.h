#ifndef SPHERE_GENERATOR_H
#define SPHERE_GENERATOR_H

#include <string>
#include <vector>
#include <memory>

#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include "vector_3d.h"

/**
 * Class for generating triangulated spherical meshes using gmsh
 * and converting them to VTK PolyData format
 */
class sphere_generator {
private:
    double radius;
    Vector3D center;
    double mesh_size;
    
    // Internal helper methods
    vtkSmartPointer<vtkPolyData> generate_sphere_with_gmsh();

public:
    /**
     * Constructor
     * @param radius Sphere radius (default: 1.0)
     * @param center Sphere center (default: origin)
     * @param mesh_size Characteristic mesh size (default: 0.1)
     */
    sphere_generator(double radius = 1.0, 
                    const Vector3D& center = Vector3D(0.0, 0.0, 0.0),
                    double mesh_size = 0.1);

    /**
     * Destructor
     */
    ~sphere_generator();

    // Getters and setters
    void set_radius(double r);
    double get_radius() const;
    
    void set_center(const Vector3D& c);
    Vector3D get_center() const;
    
    void set_mesh_size(double size);
    double get_mesh_size() const;

    /**
     * Generate triangulated spherical mesh
     * @return VTK PolyData containing the triangulated sphere
     * @throws std::runtime_error if mesh generation fails
     */
    vtkSmartPointer<vtkPolyData> generate_mesh();

    /**
     * Generate mesh and save to VTK file
     * @param output_filepath Path to save the VTK file
     * @return VTK PolyData containing the triangulated sphere
     * @throws std::runtime_error if mesh generation or file writing fails
     */
    vtkSmartPointer<vtkPolyData> generate_and_save_mesh(const std::string& output_filepath);

    /**
     * Get approximate number of triangles for current mesh size
     * @return Estimated number of triangles
     */
    size_t estimate_triangle_count() const;
};

#endif