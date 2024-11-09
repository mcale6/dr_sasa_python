#pragma once
#include "common.hpp"

class DecoupledSASA {
public:
    DecoupledSASA(float probe_radius = DEFAULT_PROBE_RADIUS, int compute_mode = 0);
    
    // Calculate from PDB file
    py::dict calculate(const std::string& pdb_file,
                      std::vector<std::vector<std::string>>& chains,
                      bool include_matrix = true);

    // Calculate from atoms
    py::dict calculate_from_atoms(std::vector<atom_struct> atoms,  // Pass by value intentionally
                                std::vector<std::vector<std::string>>& chains,
                                bool include_matrix = true);

    std::string print(std::vector<atom_struct>& atoms, const std::string& fname); // Remove const

private:
    float probe_radius_;
    int cl_mode_;
    VDWcontainer vdw_radii_;
};
