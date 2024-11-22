#include "decoupled_sasa.hpp"
#include "utils.hpp"

py::dict DecoupledSASA::calculate(const std::string& pdb_file,
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

py::dict DecoupledSASA::calculate_from_atoms(std::vector<atom_struct> atoms,
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
    
    int Imode = 2;  // Default to molecular contacts
    // Only check protein chains if no chains specified
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
            Imode = 3;
            std::vector<std::string> chain_vec(proteinChains.begin(), proteinChains.end());
            chains = {std::vector<std::string>{chain_vec[0]}, 
                     std::vector<std::string>{chain_vec[1]}};
        }
    } else {
        Imode = 3;  // Chain contacts for specified chains
    }

    if (!chains.empty()) {
        ChainSelector(chains, atoms);
    }

    SolveInteractions(atoms, Imode);
    DecoupledSolver(atoms, vdw_radii_.Points);
    CalculateDNA_ProtInteractions(atoms, 1);
    calculate_contact_areas_from_overlaps(atoms);

    py::dict results = create_analysis_results(atoms, include_matrix);
    //py::dict results; does not have sasa values, SASA = AREA - EXT0  // Total area minus buried area

    if (include_matrix) {
        results["inter_bsa_matrix"] = generate_inter_bsa_matrices(atoms);
        results["intra_bsa_matrix"] = generate_intra_bsa_matrices(atoms);
    }

    if (print_output) {
        int atmasa_bsa = 1;
        std::stringstream output;
        output << output_name;
        
        std::string overlap_file = output.str() + ".overlaps";
        PrintDNA_ProtResults(atoms, overlap_file);

        std::string dsasa_file = output.str() + ".dsasa.pdb";
        PrintDSASAResults(atoms, dsasa_file);
        
        std::string datmasa_file = output.str() + ".datmasa";
        PrintSplitAsaAtom(atoms, datmasa_file, atmasa_bsa);
        
        if(include_matrix) {
            std::string atom_matrix = output.str() + ".by_atom.tsv";
            PrintDNA_ProtResultsByAtomMatrix(atoms, atom_matrix, 1);
            
            std::string matrix_file = output.str() + ".matrix";
            Print_MatrixInsideAtom(atoms, matrix_file, 1);
        }
        
        results["printed_output"] = output.str();
    }
    return results;
}