// generic_sasa.hpp
#pragma once
#include "common.hpp"

class GenericSASA {
public:
    GenericSASA(float probe_radius = DEFAULT_PROBE_RADIUS, int compute_mode = 0);
    
    py::dict calculate(const std::string& pdb_file, 
                      std::vector<std::vector<std::string>>& chains,
                      bool include_matrix = true);

    py::dict calculate_from_atoms(std::vector<atom_struct> atoms,
                                std::vector<std::vector<std::string>>& chains,
                                bool include_matrix = true);

    std::string print(std::vector<atom_struct>& atoms, const std::string& fname);

private:
    float probe_radius_;
    int cl_mode_;
    VDWcontainer vdw_radii_;
};
