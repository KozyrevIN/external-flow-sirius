#ifndef VISUALIZER_H
#define VISUALIZER_H

#include "vector_3d.h"
#include <string>
#include <vtkDoubleArray.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>


enum class NormType {
    L1,  // L1 norm (Manhattan norm)
    L2,  // L2 norm (Euclidean norm)
    LINF // L∞ norm (Maximum norm)
};

class GradientVisualizer {
  private:
    vtkSmartPointer<vtkPolyData> mesh;
    std::string true_grad_array_name;
    std::string computed_grad_array_name;

    // Helper functions
    Vector3D getVectorFromArray(vtkDoubleArray *array, vtkIdType cellId);
    void setVectorToArray(vtkDoubleArray *array, vtkIdType cellId,
                          const Vector3D &vec);
    double computeVectorNorm(
        const Vector3D &vec,
        NormType normType); // Always uses Euclidean norm for individual vectors

  public:
    // Constructor
    GradientVisualizer(vtkSmartPointer<vtkPolyData> mesh,
                       const std::string &trueGradArrayName = "f true grad",
                       const std::string &computedGradArrayName = "Grad");

    // Main functionality methods
    void computeGradientDifference(
        const std::string &outputArrayName = "GradientDifference");

    // Computes norm over collection of per-cell Euclidean error norms
    // normType specifies how to aggregate the individual cell errors:
    // L1: sum of all cell errors, L2: sqrt(sum of squares), L∞: maximum cell
    // error
    double computeErrorNorm(NormType normType);

    // Computes Euclidean norm of error vector per cell and stores in array
    void
    computeErrorNormPerCell(NormType normType,
                            const std::string &outputArrayName = "ErrorNorm");

    // VTK visualization output methods
    void saveVisualization(const std::string &filename);

    // Utility methods
    void printErrorStatistics();
    void setTrueGradientArrayName(const std::string &name);
    void setComputedGradientArrayName(const std::string &name);
};

#endif