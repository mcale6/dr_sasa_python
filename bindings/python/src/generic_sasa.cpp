#include "generic_sasa.hpp"
#include "utils.hpp"

GenericSASA::GenericSASA(float probe_radius, int compute_mode) 
    : probe_radius_(probe_radius), cl_mode_(compute_mode) {s
    vdw_radii_.GenPoints();
}

py::dict GenericSASA::calculate(const std::string& pdb_file,
                               std::vector<std::vector<std::string>>& chains,
                               bool include_matrix,
                               bool print_output,
                               const std::string& output_name) {
    std::stringstream buffer;
    auto old_buf = std::cerr.rdbuf(buffer.rdbuf());
    auto atoms = PDBparser(pdb_file, "", true);
    if (atoms.empty()) {
        throw std::runtime_error("No atoms loaded from PDB file");
    }
    std::cerr.rdbuf(old_buf);
    
    return calculate_from_atoms(std::move(atoms), chains, include_matrix, print_output, output_name);
}

py::dict GenericSASA::calculate_from_atoms(std::vector<atom_struct> atoms,
                                         std::vector<std::vector<std::string>>& chains,
                                         bool include_matrix,
                                         bool print_output,
                                         const std::string& output_name) {
    if (atoms.empty()) {
        throw std::runtime_error("No atoms provided");
    }

    std::stringstream buffer;
    auto old_buf = std::cerr.rdbuf(buffer.rdbuf());
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

    if (!chains.empty()) {
        ChainSelector(chains, atoms);
    }

    Generic_Solver(atoms, vdw_radii_.Points, chains, Imode, cl_mode_);
    GeneratePairInteractionData(atoms);
    CalculateDNA_ProtInteractions(atoms, cl_mode_);
    
    py::dict results = create_analysis_results(atoms, include_matrix);
    
    if (print_output) {
        std::stringstream output;
        output << output_name;
        
        PrintDNA_ProtResults(atoms, output.str());
        PrintDNA_ProtResultsByAtomMatrix(atoms, output.str(), 0);
        Print_MatrixInsideAtom(atoms, output.str(), 0);
        PrintDNA_ProtResultsByAA(atoms, output.str());
        
        results["printed_output"] = output.str();
    }
    
    return results;
}
