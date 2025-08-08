#include <iostream>
#include <stdexcept>
#include "sphere_generator.h"
#include "vector_3d.h"

int main() {
    try {
        std::cout << "Sphere Generator Example\n";
        std::cout << "========================\n\n";

        // Example 1: Basic sphere generation
        std::cout << "1. Generating basic sphere (radius=1.0, mesh_size=0.2)\n";
        sphere_generator generator(1.0, Vector3D(0.0, 0.0, 0.0), 0.2);
        
        // Estimate triangle count
        size_t estimated_triangles = generator.estimate_triangle_count();
        std::cout << "   Estimated triangle count: " << estimated_triangles << std::endl;
        
        vtkSmartPointer<vtkPolyData> sphere1 = generator.generate_mesh();
        std::cout << "   Generated " << sphere1->GetNumberOfCells() << " triangles\n";
        std::cout << "   Generated " << sphere1->GetNumberOfPoints() << " points\n\n";

        // Example 2: Generate and save sphere
        std::cout << "2. Generating and saving sphere to file\n";
        vtkSmartPointer<vtkPolyData> sphere2 = generator.generate_and_save_mesh("meshes/generated_sphere.vtp");
        std::cout << "   Saved sphere to meshes/generated_sphere.vtp\n\n";

        // Example 3: Creating another sphere with different parameters
        std::cout << "3. Creating another sphere with different parameters\n";
        Vector3D center(2.0, 1.0, -0.5);
        sphere_generator small_generator(0.5, center, 0.1);
        vtkSmartPointer<vtkPolyData> sphere3 = small_generator.generate_and_save_mesh("meshes/small_sphere.vtp");
        std::cout << "   Generated small sphere at (2, 1, -0.5) with radius 0.5\n";
        std::cout << "   Generated " << sphere3->GetNumberOfCells() << " triangles\n";
        std::cout << "   Saved to meshes/small_sphere.vtp\n\n";

        // Example 4: Parameter modification
        std::cout << "4. Modifying generator parameters\n";
        generator.set_radius(1.5);
        generator.set_mesh_size(0.15);
        generator.set_center(Vector3D(1.0, 0.0, 0.0));
        
        std::cout << "   New parameters: radius=" << generator.get_radius() 
                  << ", mesh_size=" << generator.get_mesh_size() << std::endl;
        std::cout << "   Center: (" << generator.get_center().x << ", " 
                  << generator.get_center().y << ", " << generator.get_center().z << ")\n";
        
        vtkSmartPointer<vtkPolyData> sphere4 = generator.generate_mesh();
        std::cout << "   Generated " << sphere4->GetNumberOfCells() << " triangles with new parameters\n\n";

        std::cout << "All examples completed successfully!\n";
        
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        std::cerr << "\nNote: Make sure gmsh library is installed and linked properly\n";
        std::cerr << "Install gmsh development package: sudo apt install libgmsh-dev (Ubuntu/Debian)\n";
        return 1;
    }

    return 0;
}