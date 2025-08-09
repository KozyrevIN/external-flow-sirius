#include "../include/plot_generator.h"
#include "../include/geometry.h"
#include "../include/utils.h"
#include "../include/kernels.h"
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <sys/stat.h>
#include <unistd.h>
#include <map>

PlotGenerator::PlotGenerator(const std::string& plotsDirectory, 
                             const std::string& pythonScriptPath)
    : plots_dir(plotsDirectory), python_script_path(pythonScriptPath) {
    createPlotsDirectory();
    createTempDirectories();
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

void PlotGenerator::createTempDirectories() {
    // Create temp directory structure
    std::string temp_dir = "temp";
    std::string csv_dir = "temp/csv";
    std::string meshes_dir = "meshes";
    
    std::string mkdir_temp_cmd = "mkdir -p " + temp_dir;
    std::string mkdir_csv_cmd = "mkdir -p " + csv_dir;
    std::string mkdir_meshes_cmd = "mkdir -p " + meshes_dir;
    
    int result1 = std::system(mkdir_temp_cmd.c_str());
    int result2 = std::system(mkdir_csv_cmd.c_str());
    int result3 = std::system(mkdir_meshes_cmd.c_str());
    
    if (result1 == 0 && result2 == 0 && result3 == 0) {
        std::cout << "Created temporary directories: temp/csv and meshes" << std::endl;
    } else {
        std::cerr << "Warning: Some directories may not have been created properly" << std::endl;
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
    const std::string& function_name,
    const std::string& kernel_name,
    const std::string& trueGradArrayName,
    const std::string& computedGradArrayName) {
    
    try {
        // Step 1: Calculate h_max for the mesh
        double h_max = calculateMaxCellDiameter(mesh);
        std::cout << "Calculated h_max: " << h_max << std::endl;
        
        // Step 2: Create visualizer and generate CSV data
        GradientVisualizer visualizer(mesh, trueGradArrayName, computedGradArrayName);
        
        // Generate temporary CSV filename
        std::string csv_filename = "temp/csv/epsilon_analysis_" + std::to_string(getNextAvailablePlotNumber()) + ".csv";
        
        std::cout << "Generating epsilon error analysis..." << std::endl;
        visualizer.evaluateEpsilonError(
            function, kernel, epsilon_min, epsilon_max, 
            num_points, csv_filename, normType
        );
        
        // Step 3: Generate plot filename
        std::string plot_filename = getNextPlotFilename();
        
        // Step 4: Execute Python script with metadata including h_max
        std::string norm_name = (normType == NormType::L1) ? "L1" :
                               (normType == NormType::L2) ? "L2" : "L∞";
        
        std::ostringstream python_cmd;
        python_cmd << "python3 " << python_script_path 
                   << " " << csv_filename 
                   << " -o " << plot_filename
                   << " --function \"" << function_name << "\""
                   << " --kernel \"" << kernel_name << "\""
                   << " --norm \"" << norm_name << "\""
                   << " --h-max " << h_max;
        
        bool success = executeCommand(python_cmd.str());
        
        // Step 5: Clean up temporary CSV file
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
    double epsilon_max,
    const std::string& function_name,
    const std::string& kernel_name) {
    
    PlotGenerator generator;
    return generator.generateEpsilonErrorPlot(
        mesh, function, grad_function, kernel,
        epsilon_min, epsilon_max, 20, NormType::L2,
        function_name, kernel_name
    );
}

std::string PlotGenerator::generateKernelComparisonPlot(
    vtkSmartPointer<vtkPolyData> mesh,
    std::function<double(Vector3D)> function,
    std::function<Vector3D(Vector3D)> grad_function,
    double epsilon_min,
    double epsilon_max,
    int num_points,
    const std::string& function_name,
    NormType normType,
    const std::string& trueGradArrayName,
    const std::string& computedGradArrayName) {
    
    // Kernel dictionary mapping numbers to kernel functions and names
    std::map<int, std::pair<std::function<double(double)>, std::string>> kernels = {
        {2, {kernel_2, "kernel_2"}},
        {4, {kernel_4, "kernel_4"}},
        {6, {kernel_6, "kernel_6"}},
        {8, {kernel_8, "kernel_8"}}
    };
    
    try {
        // Calculate h_max for the mesh
        double h_max = calculateMaxCellDiameter(mesh);
        std::cout << "Calculated h_max: " << h_max << std::endl;
        
        // Create visualizer  
        GradientVisualizer visualizer(mesh, trueGradArrayName, computedGradArrayName);
        
        // Generate CSV files for each kernel
        std::vector<std::string> csv_files;
        std::map<int, double> best_epsilons;
        std::map<int, double> best_errors;
        
        std::cout << "Generating epsilon error analysis for all kernels..." << std::endl;
        
        for (const auto& kernel_pair : kernels) {
            int kernel_num = kernel_pair.first;
            auto kernel_func = kernel_pair.second.first;
            auto kernel_name = kernel_pair.second.second;
            
            std::string csv_filename = "temp/csv/kernel_" + std::to_string(kernel_num) + 
                                     "_analysis_" + std::to_string(getNextAvailablePlotNumber()) + ".csv";
            csv_files.push_back(csv_filename);
            
            std::cout << "Processing " << kernel_name << "..." << std::endl;
            visualizer.evaluateEpsilonError(
                function, kernel_func, epsilon_min, epsilon_max,
                num_points, csv_filename, normType
            );
            
            // Find best epsilon and error for this kernel
            std::ifstream file(csv_filename);
            std::string line;
            std::getline(file, line); // Skip header
            
            double best_error = std::numeric_limits<double>::max();
            double best_eps = 0.0;
            
            while (std::getline(file, line)) {
                std::istringstream iss(line);
                std::string epsilon_str, error_str;
                std::getline(iss, epsilon_str, ',');
                std::getline(iss, error_str, ',');
                
                double eps = std::stod(epsilon_str);
                double err = std::stod(error_str);
                
                if (err < best_error) {
                    best_error = err;
                    best_eps = eps;
                }
            }
            file.close();
            
            best_epsilons[kernel_num] = best_eps;
            best_errors[kernel_num] = best_error;
        }
        
        // Generate plot filename
        std::string plot_filename = getNextPlotFilename();
        
        // Create Python command for 2x2 subplot
        std::ostringstream python_cmd;
        python_cmd << "python3 plot_kernel_comparison.py";
        
        for (size_t i = 0; i < csv_files.size(); ++i) {
            python_cmd << " " << csv_files[i];
        }
        
        std::string norm_name = (normType == NormType::L1) ? "L1" :
                               (normType == NormType::L2) ? "L2" : "L∞";
        
        python_cmd << " -o " << plot_filename
                   << " --function \"" << function_name << "\""
                   << " --norm \"" << norm_name << "\""
                   << " --h-max " << h_max;
        
        bool success = executeCommand(python_cmd.str());
        
        // Clean up CSV files
        for (const std::string& csv_file : csv_files) {
            std::string rm_cmd = "rm -f " + csv_file;
            int rm_result = std::system(rm_cmd.c_str());
            if (rm_result != 0) {
                std::cerr << "Warning: Failed to remove temporary CSV: " << csv_file << std::endl;
            }
        }
        
        // Output best results to console
        std::cout << "\n=== KERNEL COMPARISON RESULTS ===" << std::endl;
        std::cout << "Function: " << function_name << std::endl;
        std::cout << "Epsilon range: [" << epsilon_min << ", " << epsilon_max << "]" << std::endl;
        std::cout << "Number of points: " << num_points << std::endl;
        std::cout << std::setfill('-') << std::setw(50) << "" << std::setfill(' ') << std::endl;
        
        double overall_best_error = std::numeric_limits<double>::max();
        int best_kernel = 2;
        
        for (const auto& kernel_pair : kernels) {
            int kernel_num = kernel_pair.first;
            auto kernel_name = kernel_pair.second.second;
            
            double eps = best_epsilons[kernel_num];
            double err = best_errors[kernel_num];
            
            std::cout << std::left << std::setw(10) << kernel_name 
                      << " | Best ε: " << std::setw(8) << std::fixed << std::setprecision(6) << eps
                      << " | " << norm_name << " Error: " << std::scientific << std::setprecision(3) << err << std::endl;
            
            if (err < overall_best_error) {
                overall_best_error = err;
                best_kernel = kernel_num;
            }
        }
        
        std::cout << std::setfill('-') << std::setw(50) << "" << std::setfill(' ') << std::endl;
        std::cout << "BEST METHOD: kernel_" << best_kernel 
                  << " with ε=" << std::fixed << std::setprecision(6) << best_epsilons[best_kernel]
                  << " (Error: " << std::scientific << std::setprecision(3) << overall_best_error << ")" << std::endl;
        std::cout << "==================================\n" << std::endl;
        
        if (success) {
            std::cout << "Kernel comparison plot generated: " << plot_filename << std::endl;
            return plot_filename;
        } else {
            throw std::runtime_error("Failed to execute Python comparison plotting script");
        }
        
    } catch (const std::exception& e) {
        std::cerr << "Error in generateKernelComparisonPlot: " << e.what() << std::endl;
        throw;
    }
}

std::string PlotGenerator::generate_and_linealize_plot(
    vtkSmartPointer<vtkPolyData> mesh,
    std::function<double(Vector3D)> function,
    std::function<Vector3D(Vector3D)> grad_function,
    std::function<double(double)> kernel,
    double eps_min,
    double eps_max,
    double eps_linear_min,
    double eps_linear_max,
    int num_points,
    const std::string& function_name,
    const std::string& kernel_name,
    NormType normType,
    const std::string& trueGradArrayName,
    const std::string& computedGradArrayName) {
    
    try {
        // Step 1: Generate epsilon error analysis data
        GradientVisualizer visualizer(mesh, trueGradArrayName, computedGradArrayName);
        
        std::string csv_filename = "temp/csv/linearize_analysis_" + std::to_string(getNextAvailablePlotNumber()) + ".csv";
        
        std::cout << "Generating epsilon error analysis for linearization..." << std::endl;
        visualizer.evaluateEpsilonError(
            function, kernel, eps_min, eps_max,
            num_points, csv_filename, normType
        );
        
        // Step 2: Generate plot filename
        std::string plot_filename = getNextPlotFilename();
        
        // Step 3: Execute Python script for linearization analysis
        std::string norm_name = (normType == NormType::L1) ? "L1" :
                               (normType == NormType::L2) ? "L2" : "L∞";
        
        std::ostringstream python_cmd;
        python_cmd << "python3 plot_linearization.py"
                   << " " << csv_filename
                   << " -o " << plot_filename
                   << " --eps-linear-min " << eps_linear_min
                   << " --eps-linear-max " << eps_linear_max
                   << " --function \"" << function_name << "\""
                   << " --kernel \"" << kernel_name << "\""
                   << " --norm \"" << norm_name << "\"";
        
        bool success = executeCommand(python_cmd.str());
        
        // Step 4: Clean up temporary CSV file
        std::string rm_cmd = "rm -f " + csv_filename;
        int rm_result = std::system(rm_cmd.c_str());
        if (rm_result != 0) {
            std::cerr << "Warning: Failed to remove temporary CSV: " << csv_filename << std::endl;
        }
        
        if (success) {
            std::cout << "Successfully generated linearization plot: " << plot_filename << std::endl;
            return plot_filename;
        } else {
            throw std::runtime_error("Failed to execute Python linearization plotting script");
        }
        
    } catch (const std::exception& e) {
        std::cerr << "Error in generate_and_linealize_plot: " << e.what() << std::endl;
        throw;
    }
}
//         std::cerr << "Error in generate_and_linealize_plot: " << e.what() << std::endl;
//         throw;
//     }
// }