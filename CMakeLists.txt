cmake_minimum_required(VERSION 3.10)
project(dr_sasa_simple)

# OpenMP configuration
find_package(OpenMP)
if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

# Basic configuration
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Optimization flags from original
if(NOT APPLE)
    add_definitions(-Wall -O3 -mtune=native)
    # Note: removed -mcpu=native as it might cause issues on some systems
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
)

# Set optimization flags for core library
if(OPENMP_FOUND)
    target_link_libraries(dr_sasa_core PUBLIC OpenMP::OpenMP_CXX)
endif()

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

# Print configuration info
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "OpenMP found: ${OPENMP_FOUND}")
message(STATUS "CXX Flags: ${CMAKE_CXX_FLAGS}")
message(STATUS "Python: ${Python_EXECUTABLE}")