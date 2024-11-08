// bindings/python/src/utils.hpp
#pragma once
#include "common.hpp"

py::dict create_analysis_results(const std::vector<atom_struct>& atoms, bool include_matrix = false);