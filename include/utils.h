#ifndef UTILS_H
#define UTILS_H

#include <string>

#include <vtkPolyData.h>
#include <vtkSmartPointer.h>

vtkSmartPointer<vtkPolyData> load_mesh(std::string filepath,
                                       bool verbose = false);

void save_mesh(vtkPolyData* polyData, std::string filepath);

#endif