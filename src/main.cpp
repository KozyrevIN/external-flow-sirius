#include <iostream>
#include <string>

#include <vtkPolyData.h>
#include <vtkSmartPointer.h>

#include "../include/utils.h"

int main(int argc, char *argv[]) {
    vtkSmartPointer<vtkPolyData> mesh = load_mesh(argv[1], true);

    return EXIT_SUCCESS;
}
