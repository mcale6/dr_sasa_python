#pragma once
#include "common.hpp"

class SimpleSASA {
public:
    SimpleSASA(float probe_radius = DEFAULT_PROBE_RADIUS, int compute_mode = 0);
    py::dict calculate(const std::string& pdb_file);

private:
    float probe_radius_;
    int cl_mode_;
    VDWcontainer vdw_radii_;
};