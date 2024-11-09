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
    std::cerr.rdbuf(old_buf);
    return calculate_from_atoms(std::move(atoms), chains, include_matrix);
}

py::dict GenericSASA::calculate_from_atoms(std::vector<atom_struct> atoms,
                                        std::vector<std::vector<std::string>>& chains,
                                        bool include_matrix) {
    if (atoms.empty()) {
        throw std::runtime_error("No atoms provided");
    }

    vdw_radii_.SetRadius(atoms, probe_radius_);

    int Imode;
    if (chains.empty()) {
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
        } else {
            Imode = 4;
        }
    } else {
        ChainSelector(chains, atoms);
        Imode = 1;
    }

    Generic_Solver(atoms, vdw_radii_.Points, chains, Imode, cl_mode_);
    GeneratePairInteractionData(atoms);
    CalculateDNA_ProtInteractions(atoms, cl_mode_);
    return create_analysis_results(atoms, include_matrix);
}

std::string GenericSASA::print(std::vector<atom_struct>& atoms, const std::string& fname) {
    std::stringstream output;
    output << fname;
    
    PrintDNA_ProtResults(atoms, output.str());
    PrintDNA_ProtResultsByAtomMatrix(atoms, output.str(), 0);
    Print_MatrixInsideAtom(atoms, output.str(), 0);
    PrintDNA_ProtResultsByAA(atoms, output.str());
    
    return output.str();
}