#include "sphere_generator.h"
#include "utils.h"
#include <iostream>
#include <cmath>
#include <stdexcept>
#include <map>

#include <gmsh.h>
#include <vtkPoints.h>
#include <vtkCellArray.h>
#include <vtkTriangle.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>

sphere_generator::sphere_generator(double radius, 
                                 const Vector3D& center,
                                 double mesh_size)
    : radius(radius), center(center), mesh_size(mesh_size) {
    
    if (radius <= 0.0) {
        throw std::invalid_argument("Sphere radius must be positive");
    }
    if (mesh_size <= 0.0) {
        throw std::invalid_argument("Mesh size must be positive");
    }
}

sphere_generator::~sphere_generator() {
    // Destructor - gmsh handles its own cleanup
}

void sphere_generator::set_radius(double r) {
    if (r <= 0.0) {
        throw std::invalid_argument("Sphere radius must be positive");
    }
    radius = r;
}

double sphere_generator::get_radius() const {
    return radius;
}

void sphere_generator::set_center(const Vector3D& c) {
    center = c;
}

Vector3D sphere_generator::get_center() const {
    return center;
}

void sphere_generator::set_mesh_size(double size) {
    if (size <= 0.0) {
        throw std::invalid_argument("Mesh size must be positive");
    }
    mesh_size = size;
}

double sphere_generator::get_mesh_size() const {
    return mesh_size;
}

vtkSmartPointer<vtkPolyData> sphere_generator::generate_sphere_with_gmsh() {
    // Initialize gmsh
    gmsh::initialize();
    gmsh::option::setNumber("General.Terminal", 0); // Disable terminal output
    
    try {
        // Create new model
        gmsh::model::add("sphere");
        
        // Create sphere using gmsh API
        int sphere_tag = gmsh::model::occ::addSphere(center.x, center.y, center.z, radius);
        
        // Synchronize CAD representation with gmsh model
        gmsh::model::occ::synchronize();
        
        // Set mesh size
        gmsh::vectorpair entities;
        gmsh::model::getEntities(entities, 0);
        gmsh::model::mesh::setSize(entities, mesh_size);
        
        // Generate 2D mesh (surface mesh)
        gmsh::model::mesh::generate(2);
        
        // Get mesh data from gmsh
        std::vector<std::size_t> node_tags;
        std::vector<double> node_coords;
        std::vector<double> parametric_coords;
        
        gmsh::model::mesh::getNodes(node_tags, node_coords, parametric_coords);
        
        std::vector<int> element_types;
        std::vector<std::vector<std::size_t>> element_tags;
        std::vector<std::vector<std::size_t>> element_node_tags;
        
        gmsh::model::mesh::getElements(element_types, element_tags, element_node_tags, 2); // 2D elements
        
        // Create VTK data structures
        vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
        vtkSmartPointer<vtkCellArray> triangles = vtkSmartPointer<vtkCellArray>::New();
        
        // Add points to VTK (gmsh uses 1-based indexing)
        std::map<std::size_t, vtkIdType> node_map;
        for (size_t i = 0; i < node_tags.size(); ++i) {
            vtkIdType vtk_id = points->InsertNextPoint(
                node_coords[3 * i], 
                node_coords[3 * i + 1], 
                node_coords[3 * i + 2]
            );
            node_map[node_tags[i]] = vtk_id;
        }
        
        // Add triangles to VTK
        for (size_t i = 0; i < element_types.size(); ++i) {
            if (element_types[i] == 2) { // Triangle element type in gmsh
                const std::vector<std::size_t>& nodes = element_node_tags[i];
                
                for (size_t j = 0; j < nodes.size(); j += 3) {
                    vtkSmartPointer<vtkTriangle> triangle = vtkSmartPointer<vtkTriangle>::New();
                    triangle->GetPointIds()->SetId(0, node_map[nodes[j]]);
                    triangle->GetPointIds()->SetId(1, node_map[nodes[j + 1]]);
                    triangle->GetPointIds()->SetId(2, node_map[nodes[j + 2]]);
                    triangles->InsertNextCell(triangle);
                }
            }
        }
        
        // Create VTK PolyData
        vtkSmartPointer<vtkPolyData> polyData = vtkSmartPointer<vtkPolyData>::New();
        polyData->SetPoints(points);
        polyData->SetPolys(triangles);
        
        // Clean up gmsh
        gmsh::model::remove();
        gmsh::finalize();
        
        return polyData;
        
    } catch (const std::exception& e) {
        // Clean up gmsh in case of error
        gmsh::finalize();
        throw std::runtime_error("gmsh mesh generation failed: " + std::string(e.what()));
    }
}

vtkSmartPointer<vtkPolyData> sphere_generator::generate_mesh() {
    return generate_sphere_with_gmsh();
}

vtkSmartPointer<vtkPolyData> sphere_generator::generate_and_save_mesh(const std::string& output_filepath) {
    vtkSmartPointer<vtkPolyData> mesh = generate_mesh();
    
    // Save mesh to file using existing utility function
    write_mesh(mesh, output_filepath);
    
    return mesh;
}

size_t sphere_generator::estimate_triangle_count() const {
    // Rough estimation based on sphere surface area and triangle size
    double sphere_area = 4.0 * M_PI * radius * radius;
    double triangle_area = std::sqrt(3.0) / 4.0 * mesh_size * mesh_size;
    return static_cast<size_t>(sphere_area / triangle_area);
}