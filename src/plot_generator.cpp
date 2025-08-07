#include "../include/plot_generator.h"
#include "../include/geometry.h"
#include "../include/utils.h"
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <sys/stat.h>
#include <unistd.h>

PlotGenerator::PlotGenerator(const std::string& plotsDirectory, 
                             const std::string& pythonScriptPath)
    : plots_dir(plotsDirectory), python_script_path(pythonScriptPath) {
    createPlotsDirectory();
}

void PlotGenerator::createPlotsDirectory() {
    // Use mkdir system call to create directory
    std::string mkdir_cmd = "mkdir -p " + plots_dir;
    int result = std::system(mkdir_cmd.c_str());
    if (result == 0) {
        std::cout << "Created plots directory: " << plots_dir << std::endl;
    } else {
        std::cerr << "Warning: Failed to create plots directory: " << plots_dir << std::endl;
    }
}

int PlotGenerator::getNextAvailablePlotNumber() {
    int i = 1;
    while (true) {
        std::string filename = plots_dir + "/plot_" + std::to_string(i) + ".png";
        // Check if file exists using access()
        if (access(filename.c_str(), F_OK) != 0) {
            return i;  // File doesn't exist
        }
        i++;
    }
}

std::string PlotGenerator::getNextPlotFilename() {
    int plot_number = getNextAvailablePlotNumber();
    return plots_dir + "/plot_" + std::to_string(plot_number) + ".png";
}

bool PlotGenerator::executeCommand(const std::string& command) {
    std::cout << "Executing: " << command << std::endl;
    int result = std::system(command.c_str());
    if (result == 0) {
        std::cout << "Command executed successfully." << std::endl;
        return true;
    } else {
        std::cerr << "Command failed with exit code: " << result << std::endl;
        return false;
    }
}

std::string PlotGenerator::generateEpsilonErrorPlot(
    vtkSmartPointer<vtkPolyData> mesh,
    std::function<double(Vector3D)> function,
    std::function<Vector3D(Vector3D)> grad_function,
    std::function<double(double)> kernel,
    double epsilon_min,
    double epsilon_max,
    int num_points,
    NormType normType,
    const std::string& trueGradArrayName,
    const std::string& computedGradArrayName) {
    
    try {
        // Step 1: Create visualizer and generate CSV data
        GradientVisualizer visualizer(mesh, trueGradArrayName, computedGradArrayName);
        
        // Generate temporary CSV filename
        std::string csv_filename = "temp_epsilon_analysis_" + std::to_string(getNextAvailablePlotNumber()) + ".csv";
        
        std::cout << "Generating epsilon error analysis..." << std::endl;
        visualizer.evaluateEpsilonError(
            function, kernel, epsilon_min, epsilon_max, 
            num_points, csv_filename, normType
        );
        
        // Step 2: Generate plot filename
        std::string plot_filename = getNextPlotFilename();
        
        // Step 3: Execute Python script
        std::ostringstream python_cmd;
        python_cmd << "python3 " << python_script_path 
                   << " " << csv_filename 
                   << " -o " << plot_filename;
        
        bool success = executeCommand(python_cmd.str());
        
        // Step 4: Clean up temporary CSV file
        std::string rm_cmd = "rm -f " + csv_filename;
        int rm_result = std::system(rm_cmd.c_str());
        if (rm_result == 0) {
            std::cout << "Cleaned up temporary CSV file: " << csv_filename << std::endl;
        } else {
            std::cerr << "Warning: Failed to remove temporary CSV: " << csv_filename << std::endl;
        }
        
        if (success) {
            std::cout << "Successfully generated plot: " << plot_filename << std::endl;
            return plot_filename;
        } else {
            throw std::runtime_error("Failed to execute Python plotting script");
        }
        
    } catch (const std::exception& e) {
        std::cerr << "Error in generateEpsilonErrorPlot: " << e.what() << std::endl;
        throw;
    }
}

void PlotGenerator::setPlotsDirectory(const std::string& directory) {
    plots_dir = directory;
    createPlotsDirectory();
}

void PlotGenerator::setPythonScriptPath(const std::string& path) {
    python_script_path = path;
}

std::string PlotGenerator::getPlotsDirectory() const {
    return plots_dir;
}

std::string PlotGenerator::quickPlot(
    vtkSmartPointer<vtkPolyData> mesh,
    std::function<double(Vector3D)> function,
    std::function<Vector3D(Vector3D)> grad_function,
    std::function<double(double)> kernel,
    double epsilon_min,
    double epsilon_max) {
    
    PlotGenerator generator;
    return generator.generateEpsilonErrorPlot(
        mesh, function, grad_function, kernel,
        epsilon_min, epsilon_max
    );
}