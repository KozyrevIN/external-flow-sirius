cmake_minimum_required(VERSION 3.5)

project(external-flow)

# Searching for dependencies
find_package(VTK
COMPONENTS
  CommonCore
  CommonDataModel
  IOXML
  IOLegacy
  IOCore
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

# Add sanitizers for debug builds
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -fno-omit-frame-pointer -fsanitize=address")
set(CMAKE_LINKER_FLAGS_DEBUG "${CMAKE_LINKER_FLAGS_DEBUG} -fno-omit-frame-pointer -fsanitize=address")

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