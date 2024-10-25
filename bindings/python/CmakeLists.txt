cmake_minimum_required(VERSION 3.10)
project(dr_sasa_python_bindings)

# Set visibility flags
set(CMAKE_CXX_VISIBILITY_PRESET hidden)
set(CMAKE_VISIBILITY_INLINES_HIDDEN YES)

# macOS specific settings
if(APPLE)
    set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -undefined dynamic_lookup")
endif()

# Python Configuration
find_package(Python COMPONENTS Interpreter Development REQUIRED)
find_package(pybind11 CONFIG REQUIRED)

# Create Python module
pybind11_add_module(dr_sasa_py MODULE
    bindings.cpp
)

# Include directories
target_include_directories(dr_sasa_py PRIVATE 
    ${CMAKE_SOURCE_DIR}/dr_sasa_n/src
    ${CMAKE_SOURCE_DIR}/dr_sasa_n/include
)

# Link with core library and dependencies
target_link_libraries(dr_sasa_py PRIVATE 
    dr_sasa_core
    OpenMP::OpenMP_CXX  # Add explicit OpenMP linking
)

if(USE_CUDA)
    target_link_libraries(dr_sasa_py PRIVATE ${CUDA_LIBRARIES})
endif()

if(APPLE)
    if(DEFINED LIBOMP_PREFIX)
        set_target_properties(dr_sasa_py PROPERTIES
            INSTALL_RPATH "${LIBOMP_PREFIX}/lib"
            BUILD_WITH_INSTALL_RPATH TRUE
        )
        
        target_link_libraries(dr_sasa_py PRIVATE
            ${LIBOMP_PREFIX}/lib/libomp.dylib
        )
        
        set_target_properties(dr_sasa_py PROPERTIES
            LINK_FLAGS "-Wl,-rpath,${LIBOMP_PREFIX}/lib -undefined dynamic_lookup"
        )
    endif()
endif()

# Force symbol visibility
set_target_properties(dr_sasa_py PROPERTIES 
    CXX_VISIBILITY_PRESET hidden
    VISIBILITY_INLINES_HIDDEN YES
)

# Installation configuration
if(SKBUILD)
    install(TARGETS dr_sasa_py DESTINATION .)
else()
    install(TARGETS dr_sasa_py 
            LIBRARY DESTINATION ${Python_SITEARCH}
            RUNTIME DESTINATION ${Python_SITEARCH})
endif()

# Install data files
install(FILES 
    ${CMAKE_SOURCE_DIR}/dr_sasa_n/doc/vdw.radii
    ${CMAKE_SOURCE_DIR}/dr_sasa_n/utils/points/thomson15092.xyz
    DESTINATION ${Python_SITEARCH}/dr_sasa_data
)

# Add OpenCL kernel files if they exist
if(EXISTS ${CMAKE_SOURCE_DIR}/dr_sasa_n/src/kernels)
    install(DIRECTORY ${CMAKE_SOURCE_DIR}/dr_sasa_n/src/kernels
            DESTINATION ${Python_SITEARCH}/dr_sasa_data)
endif()

message(STATUS "Current source dir: ${CMAKE_CURRENT_SOURCE_DIR}")
message(STATUS "Project source dir: ${CMAKE_SOURCE_DIR}")
message(STATUS "Python_SITEARCH: ${Python_SITEARCH}")