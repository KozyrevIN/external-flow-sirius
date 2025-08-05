cmake_minimum_required(VERSION 3.5)

project(external-flow)

# Searching for dependencies
find_package(VTK
COMPONENTS
    CommonCore
)

# Setting executable
set(SOURCES 
    src/main.cpp
)

#set(CMAKE_CXX_FLAGS "-O2")
add_executable(main ${SOURCES})

# Linking dependencies
target_link_libraries(main
  PRIVATE
   ${VTK_LIBRARIES}
)

vtk_module_autoinit(
  TARGETS main
  MODULES ${VTK_LIBRARIES}
)