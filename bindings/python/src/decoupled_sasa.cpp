#include "decoupled_sasa.hpp"
#include "utils.hpp"

DecoupledSASA::DecoupledSASA(float probe_radius, int compute_mode) 
    : probe_radius_(probe_radius), cl_mode_(compute_mode) {
    vdw_radii_.GenPoints();
}

py::dict DecoupledSASA::calculate(const std::string& pdb_file,
                                    std::vector<std::vector<std::string>>& chains,
                                    bool include_matrix) {
                                            std::stringstream buffer;
    auto old_buf = std::cerr.rdbuf(buffer.rdbuf());
    auto atoms = PDBparser(pdb_file, "", true);
    if (atoms.empty()) {
        throw std::runtime_error("No atoms loaded from PDB file");
    }
    std::cerr.rdbuf(old_buf);
    return calculate_from_atoms(std::move(atoms), chains, include_matrix);
}

py::dict DecoupledSASA::calculate_from_atoms(std::vector<atom_struct> atoms,
                                            std::vector<std::vector<std::string>>& chains,
                                            bool include_matrix) {
    if (atoms.empty()) {
        throw std::runtime_error("No atoms provided");
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

std::string DecoupledSASA::print(std::vector<atom_struct>& atoms, const std::string& fname) {
    std::stringstream output;
    output << fname;
    
    PrintDSASAResults(atoms, output.str());
    PrintDNA_ProtResultsTable(atoms, output.str());
    PrintDNA_ProtResults(atoms, output.str());
    Print_MatrixInsideAtom(atoms, output.str(), 0);
    PrintDNA_ProtResultsByAA(atoms, output.str());
    
    return output.str();
}