// bindings/python/src/decoupled_sasa.cpp
#include "decoupled_sasa.hpp"
#include "utils.hpp"

DecoupledSASA::DecoupledSASA(float probe_radius, int compute_mode) 
    : probe_radius_(probe_radius), cl_mode_(compute_mode) {
    vdw_radii_.GenPoints();
}

py::dict DecoupledSASA::calculate(const std::string& pdb_file,
                                 std::vector<std::vector<std::string>>& chains,
                                 bool include_matrix) {
    auto atoms = PDBparser(pdb_file, "", true);
    if (atoms.empty()) {
        throw std::runtime_error("No atoms loaded from PDB file");
    }

    vdw_radii_.SetRadius(atoms, probe_radius_);
    
    int Imode = (chains.size() <= 1) ? 2 : 3;

    if (!chains.empty()) {
        ChainSelector(chains, atoms);
    }

    SolveInteractions(atoms, Imode);
    DecoupledSolver(atoms, vdw_radii_.Points);
    
    return create_analysis_results(atoms, include_matrix);
}