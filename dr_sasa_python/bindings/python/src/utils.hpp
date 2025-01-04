#pragma once
#include "common.hpp"
#include "constants.hpp"

namespace py = pybind11;

// ASA component structure to hold all classifications
struct ASAComponents {
    float bb_asa = 0.0f;
    float sc_asa = 0.0f;
    float majorgroove_asa = 0.0f;
    float minorgroove_asa = 0.0f;
    float nogroove_asa = 0.0f;
    float polar_asa = 0.0f;
    float polar_bb_asa = 0.0f;
    float polar_sc_asa = 0.0f;
    float polar_majorgroove_asa = 0.0f;
    float polar_minorgroove_asa = 0.0f;
    float polar_nogroove_asa = 0.0f;
    float hyd_asa = 0.0f;
    float hyd_bb_asa = 0.0f;
    float hyd_sc_asa = 0.0f;
    float hyd_majorgroove_asa = 0.0f;
    float hyd_minorgroove_asa = 0.0f;
    float hyd_nogroove_asa = 0.0f;
    float lig_asa = 0.0f;
    float lig_polar_asa = 0.0f;
    float lig_hyd_asa = 0.0f;
};

// Other utility functions
template<typename T>
inline py::array_t<T> to_numpy(const std::vector<T>& vec) {
    return py::array_t<T>(vec.size(), vec.data());
}

template<typename T>
inline std::vector<T> from_numpy(const py::array_t<T>& arr) {
    py::buffer_info buf = arr.request();
    T* ptr = static_cast<T*>(buf.ptr);
    return std::vector<T>(ptr, ptr + buf.size);
}

// Classification helper functions
inline bool is_backbone(const std::string& atom_name) {
    return std::find(BACKBONE_ATOMS.begin(), BACKBONE_ATOMS.end(), atom_name) != BACKBONE_ATOMS.end();
}

inline bool is_polar_backbone(const std::string& atom_name) {
    return std::find(POLAR_BACKBONE.begin(), POLAR_BACKBONE.end(), atom_name) != POLAR_BACKBONE.end();
}

inline bool is_hyd_backbone(const std::string& atom_name) {
    return std::find(HYD_BACKBONE.begin(), HYD_BACKBONE.end(), atom_name) != HYD_BACKBONE.end();
}

inline bool is_polar_sidechain(const std::string& atom_name) {
    return std::find(POLAR_SIDECHAIN.begin(), POLAR_SIDECHAIN.end(), atom_name) != POLAR_SIDECHAIN.end();
}

inline bool is_hyd_sidechain(const std::string& atom_name) {
    return std::find(HYD_SIDECHAIN.begin(), HYD_SIDECHAIN.end(), atom_name) != HYD_SIDECHAIN.end();
}

inline bool is_in_major_groove(const std::string& resname, const std::string& atom_name) {
    auto it = MAJOR_GROOVE.find(resname);
    return it != MAJOR_GROOVE.end() && 
           std::find(it->second.begin(), it->second.end(), atom_name) != it->second.end();
}

inline bool is_in_minor_groove(const std::string& resname, const std::string& atom_name) {
    auto it = MINOR_GROOVE.find(resname);
    return it != MINOR_GROOVE.end() && 
           std::find(it->second.begin(), it->second.end(), atom_name) != it->second.end();
}

// Classification helper function
ASAComponents classify_atom_asa(const atom_struct& atom, float base_asa);
// Main analysis functions
py::dict create_analysis_results(const std::vector<atom_struct>& atoms, bool include_matrix);
py::dict generate_inter_bsa_matrices(std::vector<atom_struct>& atoms);
py::dict generate_intra_bsa_matrices(std::vector<atom_struct>& atoms);
void calculate_contact_areas_from_overlaps(std::vector<atom_struct>& pdb);