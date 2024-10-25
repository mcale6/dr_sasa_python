// pybind11 includes
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#include "stdafx.h"
#include "atom_struct.h"
#include "histogram.h"
#include "PDBparser2.h"
#include "SetRadius.h"
#include "SurfaceSolverOnTheFly.h"
#include "SurfaceSolverCL.h"
#include "SolverDataProcessing.h"
#include "NB.h"
#include "NonEffective.h"
#include "SearchFunctions.h"

namespace py = pybind11;

// Function we'll need based on dr_sasa.cpp
extern vector<atom_struct> PDBparser(string pdb_filename, string type_filename, bool keep_unknown);
extern void SimpleSolverCL(vector<atom_struct>& pdb, vector<float>& points, int cl_mode);

class PySASACalculator {
public:
    PySASACalculator(const std::string& vdw_file = "vdw.radii", 
                     double probe_radius = 1.4) 
        : m_rad(vdw_file),
          m_probe_radius(probe_radius),
          m_matrix_output(true)  // Default to true as in dr_sasa.cpp
    {
        m_rad.GenPoints();
    }

    // Mode 1: Chain-specific SASA with matrix output
    py::dict calculate_chain_sasa(const std::string& pdb_file, 
                                 const std::vector<std::string>& chains,
                                 const std::string& output_prefix = "") {
        try {
            auto atoms = PDBparser(pdb_file, "", true);
            
            std::vector<std::vector<std::string>> chain_sep;
            for (const auto& chain : chains) {
                chain_sep.push_back({chain});
            }
            
            if (!chains.empty()) {
                ChainSelector(chain_sep, atoms);
            }
            
            m_rad.SetRadius(atoms, m_probe_radius);
            Generic_Solver(atoms, m_rad.Points, chain_sep, 1, 0);
            GeneratePairInteractionData(atoms);

            // Prepare output files if prefix provided
            if (!output_prefix.empty()) {
                std::string datmasa = output_prefix + ".datmasa";
                std::string overlaps = output_prefix + ".overlaps";
                std::string dsasa = output_prefix + ".dsasa.pdb";

                // Output files as in dr_sasa.cpp
                PrintDNA_ProtResults(atoms, overlaps);
                PrintDSASAResults(atoms, dsasa);
                PrintSplitAsaAtom(atoms, datmasa, 1);  // Use mode 1 for bsa

                if (m_matrix_output) {
                    // Matrix outputs
                    PrintDNA_ProtResultsByAtomMatrix(atoms, output_prefix, 0);
                    Print_MatrixInsideAtom(atoms, output_prefix, 0);
                }
            }
            
            return create_results_dict(atoms);
        } catch (const std::exception& e) {
            throw std::runtime_error("Error calculating chain SASA: " + std::string(e.what()));
        }
    }

    // Mode 2: Residue dSASA mode
    py::dict calculate_residue_sasa(const std::string& pdb_file,
                                  const std::vector<std::string>& chains,
                                  const std::string& output_prefix = "") {
        try {
            auto atoms = PDBparser(pdb_file, "", true);
            
            std::vector<std::vector<std::string>> chain_sep;
            if (!chains.empty()) {
                chain_sep.push_back(chains);
                ChainSelector(chain_sep, atoms);
            }
            
            m_rad.SetRadius(atoms, m_probe_radius);
            Generic_Solver(atoms, m_rad.Points, chain_sep, 2, 0);
            GeneratePairInteractionData(atoms);

            if (!output_prefix.empty()) {
                std::string base = output_prefix + ".internal";
                if (!chains.empty()) {
                    for (const auto& c : chains) base += c;
                } else {
                    base += "all";
                }
                base += ".by_res";

                std::string int_table = base + ".int_table";
                std::string overlaps = base + ".overlaps";

                if (m_matrix_output) {
                    Print_MatrixInsideAtom(atoms, base, 0);
                } else {
                    PrintDNA_ProtResults(atoms, overlaps);
                }
            }

            return create_results_dict(atoms);
        } catch (const std::exception& e) {
            throw std::runtime_error("Error calculating residue SASA: " + std::string(e.what()));
        }
    }

    // Mode 3: Atom dSASA mode
    py::dict calculate_atom_sasa(const std::string& pdb_file,
                                const std::vector<std::string>& chains,
                                const std::string& output_prefix = "") {
        try {
            auto atoms = PDBparser(pdb_file, "", true);
            
            std::vector<std::vector<std::string>> chain_sep;
            if (!chains.empty()) {
                chain_sep.push_back(chains);
                ChainSelector(chain_sep, atoms);
            }
            
            m_rad.SetRadius(atoms, m_probe_radius);
            Generic_Solver(atoms, m_rad.Points, chain_sep, 3, 0);
            GeneratePairInteractionData(atoms);

            if (!output_prefix.empty()) {
                std::string base = output_prefix + ".internal";
                if (!chains.empty()) {
                    for (const auto& c : chains) base += c;
                } else {
                    base += "all";
                }
                base += ".by_atom";

                if (m_matrix_output) {
                    Print_MatrixInsideAtom(atoms, base, 0);
                } else {
                    PrintDNA_ProtResults(atoms, base + ".overlaps");
                }
            }

            return create_results_dict(atoms);
        } catch (const std::exception& e) {
            throw std::runtime_error("Error calculating atom SASA: " + std::string(e.what()));
        }
    }

    // Control matrix output
    void set_matrix_output(bool enable) {
        m_matrix_output = enable;
    }

    bool get_matrix_output() const {
        return m_matrix_output;
    }

private:
    // Helper to create consistent result dictionary
    py::dict create_results_dict(const std::vector<atom_struct>& atoms) {
        py::dict results;
        std::vector<double> atom_sasa;
        std::vector<std::string> atom_names;
        std::vector<std::string> residue_names;
        std::vector<int> residue_numbers;
        std::vector<std::string> chain_ids;
        
        double total_sasa = 0.0;
        for (const auto& atom : atoms) {
            atom_sasa.push_back(atom.SASA);           
            atom_names.push_back(atom.NAME);     
            residue_names.push_back(atom.RESN);
            residue_numbers.push_back(atom.RESI);
            chain_ids.push_back(atom.CHAIN);
            total_sasa += atom.SASA;
        }
        
        results["atom_sasa"] = atom_sasa;
        results["atom_names"] = atom_names;
        results["residue_names"] = residue_names;
        results["residue_numbers"] = residue_numbers;
        results["chain_ids"] = chain_ids;
        results["total_sasa"] = total_sasa;
        
        return results;
    }

    VDWcontainer m_rad;
    double m_probe_radius;
    bool m_matrix_output;
};

PYBIND11_MODULE(dr_sasa_py, m) {
    m.doc() = "Python bindings for dr_sasa: Solvent Accessible Surface Area calculator";

    py::class_<PySASACalculator>(m, "SASACalculator")
        .def(py::init<const std::string&, double>(),
             py::arg("vdw_file") = "vdw.radii",
             py::arg("probe_radius") = 1.4)
        .def("calculate_chain_sasa", &PySASACalculator::calculate_chain_sasa,
             py::arg("pdb_file"),
             py::arg("chains"),
             py::arg("output_prefix") = "",
             "Calculate chain-specific SASA with interaction matrices")
        .def("calculate_residue_sasa", &PySASACalculator::calculate_residue_sasa,
             py::arg("pdb_file"),
             py::arg("chains"),
             py::arg("output_prefix") = "",
             "Calculate residue-level dSASA")
        .def("calculate_atom_sasa", &PySASACalculator::calculate_atom_sasa,
             py::arg("pdb_file"),
             py::arg("chains"),
             py::arg("output_prefix") = "",
             "Calculate atom-level dSASA")
        .def_property("matrix_output",
                     &PySASACalculator::get_matrix_output,
                     &PySASACalculator::set_matrix_output,
                     "Enable/disable matrix output files");
}




    DrSASA(const string& vdw_file = "vdw.radii", 
           float probe_radius = 1.4,
           ComputeBackend backend = CPU)  // Add backend parameter
        : probe_radius(probe_radius),
          vdw_radii(vdw_file),
          cl_mode(0),
          matrix_output(true),
          auto_sort(true),
          compute_mode(convert_backend(backend))  // Initialize compute mode
    {
        vdw_radii.GenPoints();
        initialize_compute_device();  // Initialize GPU if available
    }
    void set_compute_backend(ComputeBackend backend) {
        compute_mode = convert_backend(backend);
        initialize_compute_device();
    }

    ComputeBackend get_compute_backend() const {
        return convert_mode(compute_mode);
    }

    string get_compute_device_info() const {
        switch(compute_mode) {
            case GPUConfig::ComputeMode::CPU:
                return "CPU (OpenMP)";
            case GPUConfig::ComputeMode::OpenCL:
                return "GPU (OpenCL)";
            case GPUConfig::ComputeMode::CUDA:
                return "GPU (CUDA)";
            default:
                return "Unknown";
        }
    }

    // Simple SASA (Mode 0)
    py::dict calculate_sasa(const string& pdb_file, 
                          bool keep_hetatm = true,
                          const string& output_prefix = "");

    // Chain Delta SASA (Mode 1)
    py::dict calculate_delta_sasa(const string& pdb_file, 
                                const vector<vector<string>>& chains,
                                const string& output_prefix = "");

    // Residue dSASA (Mode 2)
    py::dict calculate_residue_sasa(const string& pdb_file,
                                  const vector<string>& chain,
                                  const string& output_prefix = "");

    // Atom dSASA (Mode 3)
    py::dict calculate_atom_sasa(const string& pdb_file,
                               const vector<string>& chain,
                               const string& output_prefix = "");

    // Contact Surface (Mode 4)
    py::dict calculate_contact_surface(const string& pdb_file,
                                     const vector<string>& chain = {},
                                     const string& output_prefix = "");

    // Relative ASA (Mode 100)
    py::dict calculate_relative_asa(const string& pdb_file,
                                  const string& output_prefix = "");


    py::dict calculate_residue_dsasa(const string& pdb_file,
                                   const vector<string>& chains,
                                   const string& output_prefix = "") {
        atoms = PDBparser(pdb_file, "", true);
        
        if (!chains.empty()) {
            vector<vector<string>> chain_sep = {chains};
            ChainSelector(chain_sep, atoms);
        }
        
        vdw_radii.SetRadius(atoms, probe_radius);
        Generic_Solver(atoms, vdw_radii.Points, {chains}, 2, cl_mode);
        GeneratePairInteractionData(atoms);
        
        if (!output_prefix.empty()) {
            stringstream base;
            base << output_prefix << ".internal.";
            if (!chains.empty()) {
                for (const auto& c : chains) base << c;
            } else {
                base << "all";
            }
            base << ".by_res";
            
            if (matrix_output) {
                Print_MatrixInsideAtom(atoms, base.str(), 0);
            } else {
                PrintDNA_ProtResults(atoms, base.str() + ".overlaps");
            }
        }
        
        return create_residue_results_dict();
    }
    py::dict calculate_relative_asa(const string& pdb_file,
                                    const string& output_prefix = "") {
        atoms = PDBparser(pdb_file, "", true);
        vdw_radii.SetRadius(atoms, probe_radius);
        SimpleSolverCL(atoms, vdw_radii.Points, cl_mode);
        RelativeSASA(atoms);
        
        if (!output_prefix.empty()) {
            write_output(output_prefix, auto_sort ? STANDARD : ATMASA_NO_SORT);
        }
        
        return create_relative_results_dict();
    }

    // Simple SASA calculation (Mode 0)
    py::dict DrSASA::calculate_sasa(const string& pdb_file, 
                                bool keep_hetatm,
                                const string& output_prefix) {
        atoms = PDBparser(pdb_file, "", keep_hetatm);
        vdw_radii.SetRadius(atoms, probe_radius);
        SimpleSolverCL(atoms, vdw_radii.Points, cl_mode);
        
        if (!output_prefix.empty()) {
            write_output(output_prefix, auto_sort ? STANDARD : ATMASA_NO_SORT);
        }
        
        return create_results_dict();
    }

    // Delta SASA calculation (Mode 1)
    py::dict DrSASA::calculate_delta_sasa(const string& pdb_file, 
                                        const vector<vector<string>>& chains,
                                        const string& output_prefix) {
        atoms = PDBparser(pdb_file, "", true);
        
        if (!chains.empty()) {
            ChainSelector(chains, atoms);
        }
        
        vdw_radii.SetRadius(atoms, probe_radius);
        Generic_Solver(atoms, vdw_radii.Points, chains, 1, cl_mode);
        GeneratePairInteractionData(atoms);
        
        if (!output_prefix.empty()) {
            write_delta_output(output_prefix);
        }
        
        return create_interface_results_dict();
    }

    // Residue SASA calculation (Mode 2)
    py::dict DrSASA::calculate_residue_sasa(const string& pdb_file,
                                        const vector<string>& chain,
                                        const string& output_prefix) {
        atoms = PDBparser(pdb_file, "", true);
        
        if (!chain.empty()) {
            vector<vector<string>> chain_sep = {chain};
            ChainSelector(chain_sep, atoms);
        }
        
        vdw_radii.SetRadius(atoms, probe_radius);
        SimpleSolverCL(atoms, vdw_radii.Points, cl_mode);
        
        if (!output_prefix.empty()) {
            write_residue_output(output_prefix, chain);
        }
        
        return create_residue_results_dict();
    }

    // Atom SASA calculation (Mode 3)
    py::dict DrSASA::calculate_atom_sasa(const string& pdb_file,
                                        const vector<string>& chain,
                                        const string& output_prefix) {
        atoms = PDBparser(pdb_file, "", true);
        
        if (!chain.empty()) {
            vector<vector<string>> chain_sep = {chain};
            ChainSelector(chain_sep, atoms);
        }
        
        vdw_radii.SetRadius(atoms, probe_radius);
        Generic_Solver(atoms, vdw_radii.Points, {chain}, 3, cl_mode);
        GeneratePairInteractionData(atoms);
        
        if (!output_prefix.empty()) {
            write_atom_output(output_prefix, chain);
        }
        
        return create_atom_results_dict();
    }

    // Contact surface calculation (Mode 4)
    py::dict DrSASA::calculate_contact_surface(const string& pdb_file,
                                            const vector<string>& chain,
                                            const string& output_prefix) {
        atoms = PDBparser(pdb_file, "", true);
        
        if (!chain.empty()) {
            vector<vector<string>> chain_sep = {chain};
            ChainSelector(chain_sep, atoms);
        }
        
        vdw_radii.SetRadius(atoms, probe_radius);
        Generic_Solver(atoms, vdw_radii.Points, {chain}, 4, cl_mode);
        GeneratePairInteractionData(atoms);
        
        if (!output_prefix.empty()) {
            write_contact_output(output_prefix);
        }
        
        return create_contact_results_dict();
    }

    // Configuration methods
    void set_probe_radius(float radius) { probe_radius = radius; }
    float get_probe_radius() const { return probe_radius; }
    void set_matrix_output(bool enable) { matrix_output = enable; }
    bool get_matrix_output() const { return matrix_output; }
    void set_auto_sort(bool enable) { auto_sort = enable; }
    bool get_auto_sort() const { return auto_sort; }
    void set_cl_mode(int mode) { cl_mode = mode; }
    int get_cl_mode() const { return cl_mode; }
    void set_matrix_output(bool enable) { matrix_output = enable; }

private:
    struct ResidueData {
        double sasa;
        vector<const atom_struct*> atoms;
        string chain;
        int number;
        string name;
    };

    // Helper method to group atoms by residue
    map<string, ResidueData> group_by_residue() const {
        map<string, ResidueData> residues;
        
        for (const auto& atom : atoms) {
            string res_key = atom.CHAIN + "_" + 
                           atom.RESN + "_" + 
                           std::to_string(atom.RESI);
            
            if (residues.find(res_key) == residues.end()) {
                residues[res_key] = ResidueData{0.0, {}, atom.CHAIN, 
                                              atom.RESI, atom.RESN};
            }
            
            residues[res_key].sasa += atom.SASA;
            residues[res_key].atoms.push_back(&atom);
        }
        
        return residues;
    }

    static GPUConfig::ComputeMode convert_backend(ComputeBackend backend) {
        #ifdef __APPLE__
            (void)backend; // Silence unused parameter warning
            return GPUConfig::ComputeMode::CPU;  // Always return CPU for Apple
        #else
            switch(backend) {
                case CPU:
                    return GPUConfig::ComputeMode::CPU;
                case OPENCL:
                    return GPUConfig::ComputeMode::OpenCL;
                case CUDA:
                    return GPUConfig::ComputeMode::CUDA;
                default:
                    return GPUConfig::ComputeMode::CPU;
            }
        #endif
    }

    static ComputeBackend convert_mode(GPUConfig::ComputeMode mode) {
        switch(mode) {
            case GPUConfig::ComputeMode::CPU:
                return CPU;
            case GPUConfig::ComputeMode::CUDA:
                return CUDA;
            default:
                return CPU;
        }
    }

    void initialize_compute_device() {
        #ifdef __APPLE__
            compute_mode = GPUConfig::ComputeMode::CPU;
            cl_mode = 0;
        #else
            switch(compute_mode) {
                case GPUConfig::ComputeMode::CUDA:
                    cl_mode = 1;
                    break;
                default:
                    cl_mode = 0;
                    break;
            }
        #endif
    }

    // Output methods implementation
    void write_delta_output(const string& prefix) {
        string datmasa = prefix + ".datmasa";
        string overlaps = prefix + ".overlaps";
        string dsasa = prefix + ".dsasa.pdb";

        PrintDNA_ProtResults(atoms, overlaps);
        PrintDSASAResults(atoms, dsasa);
        PrintSplitAsaAtom(atoms, datmasa, auto_sort ? 1 : 3);

        if (matrix_output) {
            PrintDNA_ProtResultsByAtomMatrix(atoms, prefix, 0);
            Print_MatrixInsideAtom(atoms, prefix, 0);
        }
    }

    void write_residue_output(const string& prefix, const vector<string>& chain) {
        stringstream base;
        base << prefix << ".internal.";
        
        if (!chain.empty()) {
            for (const auto& c : chain) base << c;
        } else {
            base << "all";
        }
        base << ".by_res";

        if (matrix_output) {
            Print_MatrixInsideAtom(atoms, base.str(), 0);
        } else {
            PrintDNA_ProtResults(atoms, base.str() + ".overlaps");
        }
    }

    void write_atom_output(const string& prefix, const vector<string>& chain) {
        stringstream base;
        base << prefix << ".internal.";
        
        if (!chain.empty()) {
            for (const auto& c : chain) base << c;
        } else {
            base << "all";
        }
        base << ".by_atom";

        if (matrix_output) {
            Print_MatrixInsideAtom(atoms, base.str(), 0);
        } else {
            PrintDNA_ProtResults(atoms, base.str() + ".overlaps");
        }
    }

    void write_contact_output(const string& prefix) {
        if (matrix_output) {
            PrintDNA_ProtResultsByAtomMatrix(atoms, prefix, 1);
            Print_MatrixInsideAtom(atoms, prefix, 1);
        }
        PrintDNA_ProtResults(atoms, prefix + ".overlaps");
    }

private:
    py::dict create_results_dict() const {
        py::dict results;
        vector<double> atom_sasa;
        vector<string> atom_names;
        vector<string> residue_names;
        vector<int> residue_numbers;
        vector<string> chain_ids;
        vector<string> elements;
        vector<float> coordinates;
        
        double total_sasa = 0.0;
        for (const auto& atom : atoms) {
            atom_sasa.push_back(atom.SASA);
            atom_names.push_back(atom.NAME);
            residue_names.push_back(atom.RESN);
            residue_numbers.push_back(atom.RESI);
            chain_ids.push_back(atom.CHAIN);
            elements.push_back(atom.ELEMENT);
            coordinates.push_back(atom.COORDS[0]);
            coordinates.push_back(atom.COORDS[1]);
            coordinates.push_back(atom.COORDS[2]);
            total_sasa += atom.SASA;
        }
        
        results["atom_sasa"] = atom_sasa;
        results["atom_names"] = atom_names;
        results["residue_names"] = residue_names;
        results["residue_numbers"] = residue_numbers;
        results["chain_ids"] = chain_ids;
        results["elements"] = elements;
        results["coordinates"] = coordinates;
        results["total_sasa"] = total_sasa;
        
        return results;
    }

    py::dict create_interface_results_dict() const {
        py::dict results;
        py::list interface_atoms;
        double total_interface_area = 0.0;

        for (const auto& atom : atoms) {
            if (atom.EXT1 > 0) {  // EXT1 stores dSASA
                py::dict atom_info;
                atom_info["id"] = atom.ID;
                atom_info["name"] = atom.NAME;
                atom_info["resname"] = atom.RESN;
                atom_info["chain"] = atom.CHAIN;
                atom_info["resid"] = atom.RESI;
                atom_info["delta_sasa"] = atom.EXT1;
                atom_info["coordinates"] = py::make_tuple(atom.COORDS[0], 
                                                        atom.COORDS[1], 
                                                        atom.COORDS[2]);
                interface_atoms.append(atom_info);
                total_interface_area += atom.EXT1;
            }
        }

        results["interface_atoms"] = interface_atoms;
        results["total_interface_area"] = total_interface_area;
        return results;
    }

    py::dict create_residue_results_dict() const {
        py::dict results;
        py::list residue_list;
        double total_sasa = 0.0;
        
        auto residues = group_by_residue();
        for (const auto& [res_key, res_data] : residues) {
            py::dict res_info;
            res_info["chain"] = res_data.chain;
            res_info["number"] = res_data.number;
            res_info["name"] = res_data.name;
            res_info["sasa"] = res_data.sasa;
            
            py::list atom_list;
            for (const auto* atom : res_data.atoms) {
                py::dict atom_info;
                atom_info["name"] = atom->NAME;
                atom_info["sasa"] = atom->SASA;
                atom_list.append(atom_info);
            }
            res_info["atoms"] = atom_list;
            
            residue_list.append(res_info);
            total_sasa += res_data.sasa;
        }
        
        results["residues"] = residue_list;
        results["total_sasa"] = total_sasa;
        return results;
    }

    py::dict create_atom_results_dict() const {
        py::dict results;
        py::list atom_list;
        double total_sasa = 0.0;

        for (const auto& atom : atoms) {
            py::dict atom_info;
            atom_info["id"] = atom.ID;
            atom_info["name"] = atom.NAME;
            atom_info["resname"] = atom.RESN;
            atom_info["chain"] = atom.CHAIN;
            atom_info["resid"] = atom.RESI;
            atom_info["sasa"] = atom.SASA;
            atom_info["coordinates"] = py::make_tuple(atom.COORDS[0], 
                                                    atom.COORDS[1], 
                                                    atom.COORDS[2]);
            
            if (!atom.INTERACTION_P.empty()) {
                py::list contacts;
                for (size_t i = 0; i < atom.INTERACTION_P.size(); ++i) {
                    contacts.append(atom.INTERACTION_P[i]);
                }
                atom_info["contacts"] = contacts;
            }
            
            atom_list.append(atom_info);
            total_sasa += atom.SASA;
        }

        results["atoms"] = atom_list;
        results["total_sasa"] = total_sasa;
        return results;
    }

    py::dict create_contact_results_dict() const {
        py::dict results;
        py::list contact_list;
        double total_contact_area = 0.0;

        for (const auto& atom : atoms) {
            if (!atom.CONTACT_AREA.empty()) {
                py::dict contact_info;
                contact_info["id"] = atom.ID;
                contact_info["name"] = atom.NAME;
                contact_info["resname"] = atom.RESN;
                contact_info["chain"] = atom.CHAIN;
                contact_info["resid"] = atom.RESI;
                
                py::dict contact_areas;
                for (const auto& [contact_id, area] : atom.CONTACT_AREA) {
                    contact_areas[py::cast(contact_id)] = area;
                    total_contact_area += area;
                }
                contact_info["contact_areas"] = contact_areas;
                
                contact_list.append(contact_info);
            }
        }

        results["contacts"] = contact_list;
        results["total_contact_area"] = total_contact_area;
        return results;
    }

    py::dict create_relative_results_dict() const {
            py::dict results = create_results_dict();  // Get basic SASA results
            
            // Add relative SASA information
            vector<double> rel_sasa;
            for (const auto& atom : atoms) {
                rel_sasa.push_back(atom.SASA / atom.AREA * 100.0);  // as percentage
            }
            results["relative_sasa"] = rel_sasa;
            
            return results;
        }
    };

// Module definition
PYBIND11_MODULE(dr_sasa_py, m) {
    m.doc() = "Python bindings for dr_sasa: High-level interface";

    // Enum definition
    py::enum_<DrSASA::ComputeBackend>(m, "ComputeBackend")
        .value("CPU", DrSASA::CPU)
        .value("OPENCL", DrSASA::OPENCL)
        .value("CUDA", DrSASA::CUDA)
        .export_values();

    // Class definition
    py::class_<DrSASA>(m, "DrSASA")
        // Constructor
        .def(py::init<const string&, float, DrSASA::ComputeBackend>(),
            py::arg("vdw_file") = "vdw.radii",
            py::arg("probe_radius") = 1.4,
            py::arg("backend") = DrSASA::CPU)
    
        // Calculation methods
        .def("calculate_sasa", &DrSASA::calculate_sasa,
            py::arg("pdb_file"),
            py::arg("keep_hetatm") = true,
            py::arg("output_prefix") = "")
        .def("calculate_delta_sasa", &DrSASA::calculate_delta_sasa,
            py::arg("pdb_file"),
            py::arg("chains"),
            py::arg("output_prefix") = "")
        .def("calculate_residue_sasa", &DrSASA::calculate_residue_sasa,
            py::arg("pdb_file"),
            py::arg("chain"),
            py::arg("output_prefix") = "")
        .def("calculate_atom_sasa", &DrSASA::calculate_atom_sasa,
            py::arg("pdb_file"),
            py::arg("chain"),
            py::arg("output_prefix") = "")
        .def("calculate_contact_surface", &DrSASA::calculate_contact_surface,
            py::arg("pdb_file"),
            py::arg("chain") = vector<string>(),
            py::arg("output_prefix") = "")
        .def("calculate_relative_asa", &DrSASA::calculate_relative_asa,
            py::arg("pdb_file"),
            py::arg("output_prefix") = "")
        .def("calculate_residue_dsasa", &DrSASA::calculate_residue_dsasa,
            py::arg("pdb_file"),
            py::arg("chains"),
            py::arg("output_prefix") = "")
        // 
        .def("set_matrix_output", &DrSASA::set_matrix_output)
        .def("set_auto_sort", &DrSASA::set_auto_sort)
        // GPU configuration methods
        .def("set_compute_backend", &DrSASA::set_compute_backend)
        .def("get_compute_backend", &DrSASA::get_compute_backend)
        .def("get_compute_device_info", &DrSASA::get_compute_device_info) 
        // Properties
        .def_property("probe_radius", &DrSASA::get_probe_radius, &DrSASA::set_probe_radius)
        .def_property("matrix_output", &DrSASA::get_matrix_output, &DrSASA::set_matrix_output)
        .def_property("auto_sort", &DrSASA::get_auto_sort, &DrSASA::set_auto_sort)
        .def_property("cl_mode", &DrSASA::get_cl_mode, &DrSASA::set_cl_mode);
}