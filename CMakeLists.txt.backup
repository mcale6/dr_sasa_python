cmake_minimum_required(VERSION 3.10)
project(dr_sasa)

# Basic configuration
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_CXX_VISIBILITY_PRESET hidden)
set(CMAKE_VISIBILITY_INLINES_HIDDEN YES)

# Build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Debug")
endif()

# Compiler flags
if(APPLE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wno-c++17-extensions")
    # Define preprocessor macros to disable OpenMP and OpenCL on Apple
    add_definitions(-DNO_OPENMP -DNO_OPENCL)
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")
endif()

# Find required packages
if(NOT APPLE)
    find_package(OpenMP REQUIRED)
    find_package(OpenCL REQUIRED)
endif()

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
        ${CMAKE_SOURCE_DIR}/dr_sasa_n/include
    PRIVATE
        ${CMAKE_SOURCE_DIR}/dr_sasa_n/src
)

# Conditional linking of OpenMP and OpenCL
if(NOT APPLE)
    target_link_libraries(dr_sasa_core
        PUBLIC
            OpenMP::OpenMP_CXX
            ${OpenCL_LIBRARIES}
    )
endif()

# Python module
pybind11_add_module(dr_sasa_py MODULE
    bindings/python/bindings.cpp
)

target_include_directories(dr_sasa_py PRIVATE 
    ${CMAKE_SOURCE_DIR}/dr_sasa_n/src
    ${CMAKE_SOURCE_DIR}/dr_sasa_n/include
)

target_link_libraries(dr_sasa_py PRIVATE 
    dr_sasa_core
)

# Set output directory for the Python module
set_target_properties(dr_sasa_py PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
)

# Print debug info
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "Python: ${Python_EXECUTABLE}")
if(NOT APPLE)
    message(STATUS "OpenMP: ${OpenMP_CXX_FOUND}")
    message(STATUS "OpenCL: ${OpenCL_FOUND}")
else()
    message(STATUS "OpenMP: disabled (Apple)")
    message(STATUS "OpenCL: disabled (Apple)")
endif()
message(STATUS "Build directory: ${CMAKE_BINARY_DIR}")
message(STATUS "Install prefix: ${CMAKE_INSTALL_PREFIX}")