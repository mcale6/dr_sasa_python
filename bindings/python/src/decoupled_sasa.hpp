#pragma once
#include "common.hpp"

class DecoupledSASA {
public:
    DecoupledSASA(float probe_radius = DEFAULT_PROBE_RADIUS, int compute_mode = 0);
    
    py::dict calculate(const std::string& pdb_file,
                      std::vector<std::vector<std::string>>& chains,
                      bool include_matrix = true,
                      bool print_output = false,
                      const std::string& output_name = "output");

    py::dict calculate_from_atoms(std::vector<atom_struct> atoms,
                                std::vector<std::vector<std::string>>& chains,
                                bool include_matrix = true,
                                bool print_output = false,
                                const std::string& output_name = "output");

private:
    float probe_radius_;
    int cl_mode_;
    VDWcontainer vdw_radii_;
};