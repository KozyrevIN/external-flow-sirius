# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This is a C++ computational fluid dynamics project focused on flow around bodies and calculating surface derivatives using VTK (Visualization Toolkit). The project implements gradient calculation on mesh surfaces using kernel methods.

## Build System and Commands

### Build and Run Commands
- **Main executable**: `./scripts/run main` - builds and runs optimized version
- **Debug mode**: `./scripts/run debug` - builds with sanitizers (AddressSanitizer + UndefinedBehaviorSanitizer) and runs with memory leak detection
- **Clean build**: `./scripts/run clean` - removes build directory
- **Full rebuild**: `./scripts/run rebuild` - cleans and rebuilds all targets

### Build Requirements
- CMake 3.5+
- VTK library with components: CommonCore, CommonDataModel, IOXML, IOLegacy, FiltersGeometry
- C++17 standard
- Use WSL for building on Windows systems

### Build Targets
- `main` - optimized production build
- `debug` - debug build with sanitizers enabled and comprehensive debugging flags

## Architecture

### Core Components

**Mesh Processing Pipeline**:
1. Load mesh from VTK file (`utils.cpp::load_mesh`)
2. Attach geometric properties (cell centers, areas) (`geometry.cpp`)
3. Apply mathematical functions and compute gradients (`functions.cpp`)
4. Calculate surface derivatives using kernel methods (`geometry.cpp::grad_calculator`)

**Key Classes and Structures**:
- `Vector3D` - 3D vector structure for geometric calculations
- `grad_calculator` - Main class for computing gradients using kernel methods
  - Takes function, epsilon parameter, and kernel function
  - Implements gradient computation on mesh cells
  - Attaches computed gradients to mesh data

**Mathematical Framework**:
- Multiple kernel functions (kernel_2, kernel_4, kernel_6, kernel_8) based on Gaussian kernels
- Test functions (f_1: x-coordinate, f_2: cosine of theta angle)
- Gradient calculation using weighted kernel methods

### File Organization
- `include/` - Header files defining interfaces
- `src/` - Implementation files
- `meshes/` - Input/output mesh data (VTK format)
- `scripts/` - Build and run utilities

### VTK Integration
- Uses VTK smart pointers for automatic memory management
- Processes UnstructuredGrid and PolyData formats
- Handles mesh I/O, geometric computations, and data attachment
- **Memory Management**: When using VTK filters, explicitly clear input data (`filter->SetInputData(nullptr)`) to prevent memory leaks

## Data Flow
1. Input: VTK mesh files (sphere_642.vtk)
2. Processing: Attach centers → areas → true gradients → function values → calculated gradients
3. Output: Enhanced mesh with computed fields (sphere_642_w_fields.vtp)

## Development Notes
- Debug builds use comprehensive sanitizers - expect detailed memory analysis
- VTK memory management requires explicit cleanup in some cases
- All mesh operations work with cell-based data structures
- Kernel epsilon parameter controls the locality of gradient calculations

## Development Workflow
- **Target Management**: всегда после добавления нового таргета добавляй его в script/run