cmake_minimum_required(VERSION 3.10)
project(dr_sasa_simple)

# Basic configuration
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Find Python and pybind11
find_package(Python COMPONENTS Interpreter Development REQUIRED)
find_package(pybind11 REQUIRED)

# Core source files - only what we need for SASA
set(CORE_SOURCES
    dr_sasa_n/src/atom_struct.cpp
    dr_sasa_n/src/PDBparser2.cpp
    dr_sasa_n/src/SetRadius.cpp
    dr_sasa_n/src/SurfaceSolverCL.cpp
    dr_sasa_n/src/stdafx.cpp
    dr_sasa_n/src/SearchFunctions.cpp
)

# Core library
add_library(dr_sasa_core STATIC ${CORE_SOURCES})
target_include_directories(dr_sasa_core
    PUBLIC 
        ${CMAKE_SOURCE_DIR}/dr_sasa_n/src
)

# Python module
pybind11_add_module(dr_sasa_py MODULE
    bindings/python/bindings.cpp
)

target_include_directories(dr_sasa_py PRIVATE 
    ${CMAKE_SOURCE_DIR}/dr_sasa_n/src
)

target_link_libraries(dr_sasa_py PRIVATE 
    dr_sasa_core
)

# Set output directory
set_target_properties(dr_sasa_py PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
)