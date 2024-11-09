#pragma once
#include "common.hpp"


class SimpleSASA {
public:
    SimpleSASA(float probe_radius = DEFAULT_PROBE_RADIUS, int compute_mode = 0);
    py::dict calculate(const std::string& pdb_file);

    py::dict calculate_from_atoms(std::vector<atom_struct> atoms); 

    std::string print(std::vector<atom_struct>& atoms, const std::string& fname);

private:
    float probe_radius_;
    int cl_mode_;
    VDWcontainer vdw_radii_;
};