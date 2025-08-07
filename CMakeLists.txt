cmake_minimum_required(VERSION 3.5)

project(external-flow)

# Export compile commands for IDEs
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

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
)

#set(CMAKE_CXX_FLAGS "-O2")


add_executable(main ${SOURCES})

# Set C++ standard and include directories for main target
target_compile_features(main PRIVATE cxx_std_17)
target_include_directories(main PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)

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

target_link_libraries(debug
  PRIVATE
   ${VTK_LIBRARIES}
)

vtk_module_autoinit(
  TARGETS main debug
  MODULES ${VTK_LIBRARIES}
)