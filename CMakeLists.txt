cmake_minimum_required(VERSION 3.5)

project(external-flow)

# Export compile commands for IDEs
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Build performance optimizations
find_program(CCACHE_PROGRAM ccache)
if(CCACHE_PROGRAM)
    set(CMAKE_CXX_COMPILER_LAUNCHER "${CCACHE_PROGRAM}")
    message(STATUS "Using ccache: ${CCACHE_PROGRAM}")
endif()

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# Searching for dependencies
find_package(VTK
COMPONENTS
  CommonCore
  CommonDataModel
  IOXML
  IOLegacy
  FiltersGeometry
)

# Setting executable
set(SOURCES 
    src/main.cpp
    src/utils.cpp
    src/functions.cpp
    src/geometry.cpp
    src/visualizer.cpp
    src/plot_generator.cpp
)


set(PLOT_GENERATOR_EXAMPLE_SOURCES 
    examples/plot_generator_example.cpp
    src/utils.cpp
    src/functions.cpp
    src/geometry.cpp
    src/visualizer.cpp
    src/plot_generator.cpp
)

set(ERROR_ANALYSIS_EXAMPLE_SOURCES 
    examples/error_analysis_example.cpp
    src/utils.cpp
    src/functions.cpp
    src/geometry.cpp
    src/visualizer.cpp
    src/plot_generator.cpp
)

# Release build optimizations for main target
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -march=native")

add_executable(main ${SOURCES})
add_executable(plot_generator_example ${PLOT_GENERATOR_EXAMPLE_SOURCES})
add_executable(error_analysis_example ${ERROR_ANALYSIS_EXAMPLE_SOURCES})


# Set C++ standard and include directories for main target
target_compile_features(main PRIVATE cxx_std_17)
target_include_directories(main PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)


# Set C++ standard and include directories for plot_generator_example target
target_compile_features(plot_generator_example PRIVATE cxx_std_17)
target_include_directories(plot_generator_example PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)

# Set C++ standard and include directories for error_analysis_example target
target_compile_features(error_analysis_example PRIVATE cxx_std_17)
target_include_directories(error_analysis_example PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)

# Create debug target with maximum compatible sanitizers
add_executable(debug ${SOURCES})

# Set C++ standard and include directories for debug target
target_compile_features(debug PRIVATE cxx_std_17)
target_include_directories(debug PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)

# Add sanitizer flags for debug target - using combined approach
target_compile_options(debug PRIVATE
    -g3                           # Maximum debug info
    -O0                           # No optimization
    -fno-omit-frame-pointer      # Keep frame pointers
    -fno-optimize-sibling-calls  # Better stack traces
    -fsanitize=address,undefined # Combined AddressSanitizer and UndefinedBehaviorSanitizer
)

target_link_options(debug PRIVATE
    -fsanitize=address,undefined
)

# Linking dependencies
target_link_libraries(main
  PRIVATE
   ${VTK_LIBRARIES}
)


target_link_libraries(plot_generator_example
  PRIVATE
   ${VTK_LIBRARIES}
)

target_link_libraries(error_analysis_example
  PRIVATE
   ${VTK_LIBRARIES}
)

target_link_libraries(debug
  PRIVATE
   ${VTK_LIBRARIES}
)

vtk_module_autoinit(
  TARGETS main debug plot_generator_example error_analysis_example
  MODULES ${VTK_LIBRARIES}
)