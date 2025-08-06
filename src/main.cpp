<<<<<<< Updated upstream
#include <iostream>
#include <string>

#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkXMLPolyDataWriter.h>

#include "../include/attach_f.h"
#include "../include/functions.h"
#include "../include/kernels.h"
#include "../include/utils.h"

const std::string mesh_in_path = "meshes/sphere_642.vtk";
const std::string mesh_out_path = "meshes/sphere_642_w_fields.vtp";

int main(int argc, char *argv[]) {
    auto mesh = load_mesh(mesh_in_path, false);

    attach_f(mesh, f_2);

    auto writer = vtkSmartPointer<vtkXMLPolyDataWriter>::New();
    writer->SetFileName(mesh_out_path.c_str());
    writer->SetInputData(mesh);
    writer->Write();

    return EXIT_SUCCESS;
}
=======
#include <iostream>
#include <string>

#include <vtkPolyData.h>
#include <vtkSmartPointer.h>

#include "../include/utils.h"
#include "../include/kernels.h"
#include "../include/geometry.h"
const std::string mesh_path = "meshes/sphere_642.vtk";
#include <vtkPoints.h>

// Функция для печати координат точки по её индексу
void print_point(vtkPolyData* polyData, vtkIdType pointId) {
    vtkPoints* points = polyData->GetPoints();
    
    double point[3];
    points->GetPoint(pointId, point);
    std::cout << "Точка " << pointId << ": (" 
              << point[0] << ", " 
              << point[1] << ", " 
              << point[2] << ")" << std::endl;
}

int main(int argc, char *argv[]) {
    vtkSmartPointer<vtkPolyData> polyData = load_mesh(mesh_path,true);
    print_point(polyData, 10);
    
    // Добавляем центры ячеек к сетке
    add_center_to_cells(polyData);
    
    // Сохраняем сетку с центрами ячеек
    save_mesh(polyData, "output_with_centers.vtk");
    
    return EXIT_SUCCESS;
}
>>>>>>> Stashed changes
