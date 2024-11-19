#include "generic_sasa.hpp"
#include "utils.hpp"

GenericSASA::GenericSASA(float probe_radius, int compute_mode) 
    : probe_radius_(probe_radius), cl_mode_(compute_mode) {
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
    CalculateDNA_ProtInteractions(atoms, cl_mode_); // we are doing this twice
    
    py::dict results = create_analysis_results(atoms, include_matrix);
    
    if (include_matrix) {
        results["inter_bsa_matrix"] = generate_inter_bsa_matrices(atoms);
        results["intra_bsa_matrix"] = generate_intra_bsa_matrices(atoms);
    }

    if (print_output) {
        int atmasa_bsa = 1;
        int atmasa_sasa = 0;
        std::stringstream base_output;
        base_output << output_name;
        // .overlaps file
        std::string overlaps_file = base_output.str() + ".overlaps";
        PrintDNA_ProtResults(atoms, overlaps_file); 
        // .dsasa.pdb file
        std::string dsasa_file = base_output.str() + ".dsasa.pdb";
        PrintDSASAResults(atoms, dsasa_file);
        // .datmasa file
        std::string datmasa_file = base_output.str() + ".datmasa";
        PrintSplitAsaAtom(atoms, datmasa_file, atmasa_bsa);
        // .asa.pdb file - main SASA results in PDB format
        std::string asa_file = base_output.str() + ".asa.pdb";
        PrintSASAResults(atoms, asa_file);
        // .atmasa file - per-type SASA analysis
        std::string atmasa_file = base_output.str() + ".atmasa";
        PrintSplitAsaAtom(atoms, atmasa_file, atmasa_sasa);
        if(include_matrix) {
            // By atom matrix files
            std::string atom_matrix_file = base_output.str() + ".by_atom.tsv";
            PrintDNA_ProtResultsByAtomMatrix(atoms, atom_matrix_file, 0);
            // Matrix inside atom files
            std::string matrix_inside_file = base_output.str() + ".matrix";
            Print_MatrixInsideAtom(atoms, matrix_inside_file, 0);
        }
        results["printed_output"] = base_output.str();
    }
    
    return results;
}