# Minimum CMake version
cmake_minimum_required(VERSION 3.10)

# Project definition
project(dr_sasa_python VERSION 0.1.0)

# Basic configuration
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Virtual Environment Detection
if(DEFINED ENV{VIRTUAL_ENV})
    message(STATUS "Virtual environment found at: $ENV{VIRTUAL_ENV}")
    set(Python_ROOT_DIR "$ENV{VIRTUAL_ENV}")
    set(Python_EXECUTABLE "$ENV{VIRTUAL_ENV}/bin/python")
else()
    message(WARNING "No virtual environment detected! It's recommended to use a virtual environment.")
endif()

# Platform-specific compiler settings
if(APPLE)
    set(CMAKE_C_COMPILER "/opt/homebrew/opt/llvm/bin/clang")
    set(CMAKE_CXX_COMPILER "/opt/homebrew/opt/llvm/bin/clang++")
else()
    # Optimization flags for non-Apple platforms
    add_definitions(-Wall -O3 -march=native -mtune=native)
endif()

# Python configuration
find_package(Python COMPONENTS Interpreter Development REQUIRED)
# Installation prefix configuration
if(SKBUILD)
    set(CMAKE_INSTALL_PREFIX ${SKBUILD_PLATLIB}/dr_sasa_python)
    set(INSTALL_DIR "bindings/python")  # Relative to dr_sasa_python package
else()
    if(DEFINED ENV{VIRTUAL_ENV})
        # Make sure to find Python before using its version
        find_package(Python COMPONENTS Interpreter Development REQUIRED)
        # Use explicit version numbers
        set(CMAKE_INSTALL_PREFIX "$ENV{VIRTUAL_ENV}/lib/python${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}/site-packages/dr_sasa_python")
        message(STATUS "Setting install prefix to virtual environment: ${CMAKE_INSTALL_PREFIX}")
    else()
        set(CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}/lib/python${Python_VERSION_MAJOR}.${Python_VERSION_MINOR}/site-packages/dr_sasa_python")
    endif()
    set(INSTALL_DIR "bindings/python")  # Keep same structure
endif()

# Dependencies
## 1. OpenMP Configuration
if(APPLE)
    # OpenMP Settings for Apple
    set(OpenMP_C_FLAGS "-Xpreprocessor -fopenmp")
    set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp")
    set(OpenMP_C_LIB_NAMES "omp")
    set(OpenMP_CXX_LIB_NAMES "omp")
    set(OpenMP_omp_LIBRARY "/opt/homebrew/opt/libomp/lib/libomp.dylib")
    include_directories("/opt/homebrew/opt/libomp/include")
endif()

find_package(OpenMP)
if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -g")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

## 2. OpenCL Configuration
find_package(OpenCL REQUIRED)
if(OpenCL_FOUND)
    message(STATUS "OpenCL found: ${OpenCL_LIBRARIES}")
    add_definitions(-DCL_TARGET_OPENCL_VERSION=120)
    add_definitions(-DUSE_GPU)
endif()

## 3. pybind11
find_package(pybind11 REQUIRED)

# Source files configuration
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

# Library targets
## 1. Core library
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

## 2. Python module
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

# Installation configuration
install(TARGETS dr_sasa_py 
        LIBRARY DESTINATION ${INSTALL_DIR}
        RUNTIME DESTINATION ${INSTALL_DIR})

# Package version configuration
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    "${CMAKE_BINARY_DIR}/dr_sasa_pythonConfigVersion.cmake"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY AnyNewerVersion
)

# Configuration summary
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "Installation prefix: ${CMAKE_INSTALL_PREFIX}")
message(STATUS "Install directory: ${INSTALL_DIR}")
message(STATUS "OpenMP found: ${OPENMP_FOUND}")
message(STATUS "OpenCL found: ${OpenCL_FOUND}")
message(STATUS "CXX Flags: ${CMAKE_CXX_FLAGS}")
message(STATUS "Python: ${Python_EXECUTABLE}")