#include "simple_sasa.hpp"
#include "utils.hpp"

SimpleSASA::SimpleSASA(float probe_radius, int compute_mode) 
    : probe_radius_(probe_radius), cl_mode_(compute_mode) {
    vdw_radii_.GenPoints();
}

py::dict SimpleSASA::calculate(const std::string& pdb_file,
                             bool print_output,
                             const std::string& output_name) {
    std::stringstream buffer;
    auto old_buf = std::cerr.rdbuf(buffer.rdbuf());
    auto atoms = PDBparser(pdb_file, "", true);
    if (atoms.empty()) {
        throw std::runtime_error("No atoms loaded from PDB file");
    }
    std::cerr.rdbuf(old_buf);
    
    return calculate_from_atoms(std::move(atoms), print_output, output_name);
}

py::dict SimpleSASA::calculate_from_atoms(std::vector<atom_struct> atoms,
                                        bool print_output,
                                        const std::string& output_name) {
    if (atoms.empty()) {
        throw std::runtime_error("No atoms provided");
    }

    std::stringstream buffer;
    auto old_buf = std::cerr.rdbuf(buffer.rdbuf());
    vdw_radii_.SetRadius(atoms, probe_radius_);
    std::cerr.rdbuf(old_buf);
    
    SolveInteractions(atoms, 0);
    SimpleSolverCL(atoms, vdw_radii_.Points, cl_mode_);
    
    py::dict results = create_analysis_results(atoms, false);
    
    if (print_output) {
        int atmasa_sasa = 0;
        std::stringstream base_output;
        base_output << output_name;

        // .asa.pdb file - main SASA results in PDB format
        std::string asa_file = base_output.str() + ".asa.pdb";
        PrintSASAResults(atoms, asa_file);
        // .atmasa file - per-type SASA analysis
        std::string atmasa_file = base_output.str() + ".atmasa";
        PrintSplitAsaAtom(atoms, atmasa_file, atmasa_sasa);
        
        results["printed_output"] = base_output.str();
    }
    
    return results;
}