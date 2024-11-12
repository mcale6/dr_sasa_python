#pragma once
#include "common.hpp"

void RelativeSASA2(std::vector<atom_struct>& atoms);

py::dict create_analysis_results(const std::vector<atom_struct>& atoms, bool include_matrix);

py::dict create_analysis_results(const std::vector<atom_struct>& atoms, bool include_matrix);

py::dict generate_interaction_matrices(std::vector<atom_struct>& atoms);

py::dict generate_intra_matrices(std::vector<atom_struct>& atoms);

namespace conversion {

// Convert 1D vector representing 2D matrix to 2D numpy array
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

// Optional convenience function for atom contacts
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

}