// bindings/python/src/generic_sasa.cpp
#include "generic_sasa.hpp"
#include "utils.hpp"

GenericSASA::GenericSASA(float probe_radius, int compute_mode) 
    : probe_radius_(probe_radius), cl_mode_(compute_mode) {
    vdw_radii_.GenPoints();
}

py::dict GenericSASA::calculate(const std::string& pdb_file, 
                               std::vector<std::vector<std::string>>& chains,
                               bool include_matrix) {
    std::stringstream buffer;
    auto old_buf = std::cerr.rdbuf(buffer.rdbuf());
    auto atoms = PDBparser(pdb_file, "", true);
    if (atoms.empty()) {
        throw std::runtime_error("No atoms loaded from PDB file");
    }
    vdw_radii_.SetRadius(atoms, probe_radius_);
    std::cerr.rdbuf(old_buf);

    int Imode;
    if (chains.size() <= 1) {
        std::set<std::string> proteinChains;
        for (const auto& atom: atoms) {
            if (atom.MOL_TYPE == "PROTEIN") {
                proteinChains.insert(atom.CHAIN);
            } else {
                proteinChains.clear();
                break;
            }
        }
        if (proteinChains.size() == 2) {
            Imode = 5;
            std::vector<std::string> chain_vec(proteinChains.begin(), proteinChains.end());
            chains = {std::vector<std::string>{chain_vec[0]}, std::vector<std::string>{chain_vec[1]}};
        } else {
            Imode = 4;
        }
    } else {
        Imode = 1;
    }

    ChainSelector(chains, atoms);
    Generic_Solver(atoms, vdw_radii_.Points, chains, Imode, cl_mode_);
    GeneratePairInteractionData(atoms);
    CalculateDNA_ProtInteractions(atoms, cl_mode_);
    return create_analysis_results(atoms, include_matrix);
}