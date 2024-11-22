#pragma once
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <stdexcept>
#include <cstdint>
#include <algorithm>
#include <cmath>
#include <filesystem>
#include <sstream>

#include "stdafx.h"
#include "atom_struct.h"
#include "PDBparser2.h"
#include "SetRadius.h"
#include "SurfaceSolverCL.h"
#include "SurfaceSolverOnTheFly.h"
#include "histogram.h"
#include "SolverDataProcessing.h"

namespace py = pybind11;
using namespace py::literals;

// Utility functions
template<typename T>
inline py::array_t<T> vector_to_numpy(const std::vector<T>& vec) {
    return py::array_t<T>(vec.size(), vec.data());
}

template<typename T>
inline std::vector<T> numpy_to_vector(const py::array_t<T>& arr) {
    py::buffer_info buf = arr.request();
    T* ptr = static_cast<T*>(buf.ptr);
    return std::vector<T>(ptr, ptr + buf.size);
}

// Convert 1D vector representing 2D matrix to 2D numpy array
template<typename T>
inline py::array_t<float> contact_matrix_to_numpy(
    const std::vector<float>& matrix,
    size_t rows,
    size_t cols
) {
    if (matrix.size() != rows * cols) {
        throw std::runtime_error("Matrix size doesn't match dimensions");
    }
    return py::array_t<float>({rows, cols}, matrix.data());
}

template<typename T>
inline py::dict atom_contacts_to_dict(const atom_struct& atom) {
    py::dict result;
    for (const auto& [id, area] : atom.CONTACT_AREA) {
        result[py::cast(id)] = py::dict(
            "area"_a=area,
            "distance"_a=atom.DISTANCES.count(id) ? atom.DISTANCES.at(id) : 0.0f
        );
    }
    return result;
}