#pragma once
#include "common.hpp"

class SimpleSASA {
public:
    // Constructor with member initialization list
    SimpleSASA(float probe_radius = DEFAULT_PROBE_RADIUS, int compute_mode = 0)
        : vdw_radii_(""),  // Initialize VDWcontainer directly
          probe_radius_(probe_radius),
          cl_mode_(compute_mode)
    {
        vdw_radii_.GenPoints();
    }
    
    py::dict calculate(const std::string& pdb_file,
                      bool print_output = false,
                      const std::string& output_name = "output");

    py::dict calculate_from_atoms(std::vector<atom_struct> atoms,
                                bool print_output = false,
                                const std::string& output_name = "output");

private:
    VDWcontainer vdw_radii_;  // Now initialized in constructor
    float probe_radius_;
    int cl_mode_;
};