cmake_minimum_required(VERSION 3.10)
project(dr_sasa_python VERSION 0.1.0)

# Set installation prefix and paths
if(SKBUILD)
    set(CMAKE_INSTALL_PREFIX ${SKBUILD_PLATLIB}/dr_sasa_python)
    set(INSTALL_DIR "bindings/python")  # Relative to dr_sasa_python package
else()
    set(CMAKE_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX}/lib/python${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}/site-packages/dr_sasa_python)
    set(INSTALL_DIR "bindings/python")  # Keep same structure
endif()

# OpenMP configuration
find_package(OpenMP)
if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -g")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

# Basic configuration
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Optimization flags
if(NOT APPLE)
    add_definitions(-Wall -O3 -march=native -mtune=native)
endif()

# OpenCL configuration
find_package(OpenCL REQUIRED)
if(OpenCL_FOUND)
    message(STATUS "OpenCL found: ${OpenCL_LIBRARIES}")
    add_definitions(-DCL_TARGET_OPENCL_VERSION=120)
    add_definitions(-DUSE_GPU)
endif()

# Find Python and pybind11
find_package(Python COMPONENTS Interpreter Development REQUIRED)
find_package(pybind11 REQUIRED)

# Core source files
set(CORE_SOURCES
    dr_sasa_n/src/atom_struct.cpp
    dr_sasa_n/src/residue_struct.cpp
    dr_sasa_n/src/PDBparser2.cpp
    dr_sasa_n/src/SetRadius.cpp
    dr_sasa_n/src/SurfaceSolverCL.cpp
    dr_sasa_n/src/NB.cpp
    dr_sasa_n/src/SolverDataProcessing.cpp
    dr_sasa_n/src/NonEffective.cpp
    dr_sasa_n/src/histogram.cpp
    dr_sasa_n/src/SearchFunctions.cpp
    dr_sasa_n/src/SurfaceSolverOnTheFly.cpp
    dr_sasa_n/src/dr_sasa.cpp
    dr_sasa_n/src/stdafx.cpp
)

# Core library
add_library(dr_sasa_core STATIC ${CORE_SOURCES})
target_include_directories(dr_sasa_core
    PUBLIC 
        ${CMAKE_SOURCE_DIR}/dr_sasa_n/src
        ${OpenCL_INCLUDE_DIRS}
)

target_link_libraries(dr_sasa_core
    PUBLIC 
        $<$<BOOL:${OPENMP_FOUND}>:OpenMP::OpenMP_CXX>
        $<$<BOOL:${OpenCL_FOUND}>:OpenCL::OpenCL>
)

# Python module (creates dr_sasa_py.so)
pybind11_add_module(dr_sasa_py MODULE
    dr_sasa_python/bindings/python/src/atom_bindings.cpp
    dr_sasa_python/bindings/python/src/sasa_calculators.cpp
    dr_sasa_python/bindings/python/src/utils.cpp
    dr_sasa_python/bindings/python/bindings.cpp
)

target_include_directories(dr_sasa_py PRIVATE 
    ${CMAKE_SOURCE_DIR}/dr_sasa_n/src
    ${OpenCL_INCLUDE_DIRS}
)

target_link_libraries(dr_sasa_py PRIVATE 
    dr_sasa_core
)

# Installation - puts dr_sasa_py.so in the right place
install(TARGETS dr_sasa_py 
        LIBRARY DESTINATION ${INSTALL_DIR}
        RUNTIME DESTINATION ${INSTALL_DIR})

# Package version file
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    "${CMAKE_BINARY_DIR}/dr_sasa_pythonConfigVersion.cmake"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY AnyNewerVersion
)

# Print configuration info
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "Installation prefix: ${CMAKE_INSTALL_PREFIX}")
message(STATUS "Install directory: ${INSTALL_DIR}")
message(STATUS "OpenMP found: ${OPENMP_FOUND}")
message(STATUS "OpenCL found: ${OpenCL_FOUND}")
message(STATUS "CXX Flags: ${CMAKE_CXX_FLAGS}")
message(STATUS "Python: ${Python_EXECUTABLE}")