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
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")
endif()

# OpenMP Configuration
if(APPLE)
    execute_process(
        COMMAND brew --prefix libomp
        OUTPUT_VARIABLE LIBOMP_PREFIX
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    if(LIBOMP_PREFIX)
        set(OpenMP_C_FLAGS "-Xpreprocessor -fopenmp -I${LIBOMP_PREFIX}/include")
        set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp -I${LIBOMP_PREFIX}/include")
        set(OpenMP_C_LIB_NAMES "omp")
        set(OpenMP_CXX_LIB_NAMES "omp")
        set(OpenMP_omp_LIBRARY "${LIBOMP_PREFIX}/lib/libomp.dylib")
        link_directories("${LIBOMP_PREFIX}/lib")
        set(CMAKE_INSTALL_RPATH "${LIBOMP_PREFIX}/lib")
        set(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE)
    else()
        message(FATAL_ERROR "OpenMP not found. Install it with 'brew install libomp'")
    endif()
endif()

# Find required packages
find_package(OpenMP REQUIRED)
find_package(OpenCL REQUIRED)
find_package(Python COMPONENTS Interpreter Development REQUIRED)
find_package(pybind11 REQUIRED)

# Core source files
set(CORE_SOURCES
    dr_sasa_n/src/stdafx.cpp
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
    dr_sasa_n/src/dr_sasa.cpp
    dr_sasa_n/src/SurfaceSolverOnTheFly.cpp
)

# Core library
add_library(dr_sasa_core STATIC ${CORE_SOURCES})
target_include_directories(dr_sasa_core
    PUBLIC 
        ${CMAKE_SOURCE_DIR}/dr_sasa_n/include
    PRIVATE
        ${CMAKE_SOURCE_DIR}/dr_sasa_n/src
)
target_link_libraries(dr_sasa_core
    PUBLIC
        OpenMP::OpenMP_CXX
        ${OpenCL_LIBRARIES}
)

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
message(STATUS "OpenMP: ${OpenMP_CXX_FOUND}")
message(STATUS "OpenCL: ${OpenCL_FOUND}")
message(STATUS "Build directory: ${CMAKE_BINARY_DIR}")
message(STATUS "Install prefix: ${CMAKE_INSTALL_PREFIX}")