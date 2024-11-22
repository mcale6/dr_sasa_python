#pragma once
#include "common.hpp"

class DecoupledSASA {
public:
    // Constructor with member initialization list
    DecoupledSASA(float probe_radius = DEFAULT_PROBE_RADIUS, int compute_mode = 0)
        : vdw_radii_(""),  // Initialize VDWcontainer directly
          probe_radius_(probe_radius),
          cl_mode_(compute_mode)
    {
        vdw_radii_.GenPoints();
    }
    
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
    VDWcontainer vdw_radii_;
    float probe_radius_;
    int cl_mode_;
};