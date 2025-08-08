
#ifndef PLOT_GENERATOR_H
#define PLOT_GENERATOR_H

#include "visualizer.h"
#include "vector_3d.h"
#include <functional>
#include <string>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>

class PlotGenerator {
private:
    std::string plots_dir;
    std::string python_script_path;
    
    // Helper functions
    void createPlotsDirectory();
    void createTempDirectories();
    std::string getNextPlotFilename();
    int getNextAvailablePlotNumber();
    bool executeCommand(const std::string& command);
    
public:
    // Constructor
    PlotGenerator(const std::string& plotsDirectory = "plots", 
                  const std::string& pythonScriptPath = "plot_epsilon_error.py");
    
    // Main functionality
    std::string generateEpsilonErrorPlot(
        vtkSmartPointer<vtkPolyData> mesh,
        std::function<double(Vector3D)> function,
        std::function<Vector3D(Vector3D)> grad_function,
        std::function<double(double)> kernel,
        double epsilon_min = 0.01,
        double epsilon_max = 1.0,
        int num_points = 20,
        NormType normType = NormType::L2,
        const std::string& function_name = "f",
        const std::string& kernel_name = "kernel",
        const std::string& trueGradArrayName = "f true grad",
        const std::string& computedGradArrayName = "Grad"
    );
    
    // Utility methods
    void setPlotsDirectory(const std::string& directory);
    void setPythonScriptPath(const std::string& path);
    std::string getPlotsDirectory() const;
    
    // Static helper for quick plot generation
    static std::string quickPlot(
        vtkSmartPointer<vtkPolyData> mesh,
        std::function<double(Vector3D)> function,
        std::function<Vector3D(Vector3D)> grad_function,
        std::function<double(double)> kernel,
        double epsilon_min = 0.01,
        double epsilon_max = 1.0,
        const std::string& function_name = "f",
        const std::string& kernel_name = "kernel"
    );
    
    // Multi-kernel comparison plot
    std::string generateKernelComparisonPlot(
        vtkSmartPointer<vtkPolyData> mesh,
        std::function<double(Vector3D)> function,
        std::function<Vector3D(Vector3D)> grad_function,
        double epsilon_min = 0.01,
        double epsilon_max = 1.0,
        int num_points = 20,
        const std::string& function_name = "f",
        NormType normType = NormType::L2,
        const std::string& trueGradArrayName = "f true grad",
        const std::string& computedGradArrayName = "Grad"
    );

    // Generate plot with linear regression analysis
    std::string generate_and_linealize_plot(
        vtkSmartPointer<vtkPolyData> mesh,
        std::function<double(Vector3D)> function,
        std::function<Vector3D(Vector3D)> grad_function,
        std::function<double(double)> kernel,
        double eps_min,
        double eps_max,
        double eps_linear_min,
        double eps_linear_max,
        int num_points = 50,
        const std::string& function_name = "f",
        const std::string& kernel_name = "kernel",
        NormType normType = NormType::L2,
        const std::string& trueGradArrayName = "f true grad",
        const std::string& computedGradArrayName = "Grad"
    );
    
};

#endif