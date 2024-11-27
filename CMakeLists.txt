cmake_minimum_required(VERSION 3.10)
project(dr_sasa_simple)

# OpenMP configuration with GPU target offloading
find_package(OpenMP)
if(OPENMP_FOUND)
    # Check for NVPTX offloading support
    find_program(NVPTX_OFFLOAD_COMPILER "nvptx-none-gcc")
    if(NVPTX_OFFLOAD_COMPILER)
        message(STATUS "NVPTX offloading compiler found at: ${NVPTX_OFFLOAD_COMPILER}")
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS} -fopenmp-targets=nvptx64-nvidia-cuda")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -fopenmp-targets=nvptx64-nvidia-cuda -g")
        add_definitions(-DUSE_OPENMP_GPU_OFFLOAD)
    else()
        message(WARNING "NVPTX offloading compiler not found. Install gcc-offload-nvptx package for GPU support.")
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -g")
    endif()
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

# Basic configuration
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Optimization flags from original
if(NOT APPLE)
    add_definitions(-Wall -O3 -march=native -mtune=native)
endif()

# OpenCL configuration
find_package(OpenCL REQUIRED)
if(OpenCL_FOUND)
    message(STATUS "OpenCL found: ${OpenCL_LIBRARIES}")
    add_definitions(-DCL_TARGET_OPENCL_VERSION=120)
    add_definitions(-DUSE_GPU)
else()
    message(WARNING "OpenCL not found, OpenCL support will be disabled")
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

# Copy OpenCL kernel files to build directory
file(GLOB CL_KERNELS "dr_sasa_n/src/*.cl")
file(COPY ${CL_KERNELS} DESTINATION ${CMAKE_BINARY_DIR}/kernels)

# Core library
add_library(dr_sasa_core STATIC ${CORE_SOURCES})
target_include_directories(dr_sasa_core
    PUBLIC 
        ${CMAKE_SOURCE_DIR}/dr_sasa_n/src
        ${OpenCL_INCLUDE_DIRS}
)

# Link OpenCL and OpenMP
target_link_libraries(dr_sasa_core
    PUBLIC 
        $<$<BOOL:${OPENMP_FOUND}>:OpenMP::OpenMP_CXX>
        $<$<BOOL:${OpenCL_FOUND}>:OpenCL::OpenCL>
)

# Python module
pybind11_add_module(dr_sasa_py MODULE
    bindings/python/src/atom_bindings.cpp
    bindings/python/src/sasa_calculators.cpp
    bindings/python/src/utils.cpp
    bindings/python/bindings.cpp
)

target_include_directories(dr_sasa_py PRIVATE 
    ${CMAKE_SOURCE_DIR}/dr_sasa_n/src
    ${OpenCL_INCLUDE_DIRS}
)

target_link_libraries(dr_sasa_py PRIVATE 
    dr_sasa_core
)

# Set output directory
set_target_properties(dr_sasa_py PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
)

# Copy OpenCL kernels to installation directory
install(DIRECTORY ${CMAKE_BINARY_DIR}/kernels/
        DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/kernels)

# Print configuration info
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")
message(STATUS "OpenMP found: ${OPENMP_FOUND}")
if(OPENMP_FOUND)
    message(STATUS "OpenMP GPU offload: ${OPENMP_HAVE_NVPTX_OFFLOAD}")
endif()
message(STATUS "OpenCL found: ${OpenCL_FOUND}")
message(STATUS "CXX Flags: ${CMAKE_CXX_FLAGS}")
message(STATUS "Python: ${Python_EXECUTABLE}")