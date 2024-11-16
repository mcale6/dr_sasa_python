#include "decoupled_sasa.hpp"
#include "utils.hpp"

DecoupledSASA::DecoupledSASA(float probe_radius, int compute_mode) 
    : probe_radius_(probe_radius), cl_mode_(compute_mode) {
    vdw_radii_.GenPoints();
}

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
    
    int Imode = (chains.size() <= 1) ? 2 : 3;

    if (!chains.empty()) {
        ChainSelector(chains, atoms);
    }

    SolveInteractions(atoms, Imode);
    DecoupledSolver(atoms, vdw_radii_.Points); // in the original implementation the results are not saved in atom struct.  
    
    py::dict results = create_analysis_results(atoms, include_matrix);
    

    if (include_matrix) {
        results["inter_bsa_matrix"] = generate_inter_bsa_matrices(atoms);
        results["intra_bsa_matrix"] = generate_intra_bsa_matrices(atoms);
    }
    
    if (print_output) {
        std::stringstream output;
        output << output_name;
        
        PrintDSASAResults(atoms, output.str());
        PrintDNA_ProtResultsTable(atoms, output.str());
        PrintDNA_ProtResults(atoms, output.str());
        Print_MatrixInsideAtom(atoms, output.str(), 0);
        PrintDNA_ProtResultsByAA(atoms, output.str());
        
        results["printed_output"] = output.str();
    }
    
    return results;
}