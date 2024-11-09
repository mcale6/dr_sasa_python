#include "simple_sasa.hpp"
#include "utils.hpp"

SimpleSASA::SimpleSASA(float probe_radius, int compute_mode) 
    : probe_radius_(probe_radius), cl_mode_(compute_mode) {
    vdw_radii_.GenPoints();
}

py::dict SimpleSASA::calculate(const std::string& pdb_file) {
    std::stringstream buffer;
    auto old_buf = std::cerr.rdbuf(buffer.rdbuf());
    auto atoms = PDBparser(pdb_file, "", true);
    if (atoms.empty()) {
        throw std::runtime_error("No atoms loaded from PDB file");
    }
    std::cerr.rdbuf(old_buf);
    
    return calculate_from_atoms(std::move(atoms));
}

py::dict SimpleSASA::calculate_from_atoms(std::vector<atom_struct> atoms) {
    // Set radius
    vdw_radii_.SetRadius(atoms, probe_radius_);
    
    // Calculate SASA
    SolveInteractions(atoms, 0);
    SimpleSolverCL(atoms, vdw_radii_.Points, cl_mode_);

    return create_analysis_results(atoms);
}

std::string SimpleSASA::print(std::vector<atom_struct>& atoms, const std::string& fname) {
    std::stringstream output;
    output << fname;
    
    PrintSASAResults(atoms, output.str());
    PrintSplitAsaAtom(atoms, output.str(), 0);
    
    return output.str();
}