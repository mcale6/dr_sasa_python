// bindings/python/src/common.hpp
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

static constexpr float DEFAULT_PROBE_RADIUS = 1.4f;

// Common exceptions
class SASAError : public std::runtime_error {
    using std::runtime_error::runtime_error;
};

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