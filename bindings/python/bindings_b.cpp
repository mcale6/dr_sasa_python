#include <string>
#include <vector>
#include <cstdint>
#include <map>
#include <memory>
#include <stdexcept>
#include <sstream>
// pybind11 headers
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
// Project headers
#include "stdafx.h"
#include "atom_struct.h"
#include "PDBparser2.h"
#include "SetRadius.h"
#include "histogram.h"
#include "SurfaceSolverOnTheFly.h"
#include "SurfaceSolverCL.h"
#include "SolverDataProcessing.h"
#include "NB.h"
#include "gpu_config.hpp"

using std::string;
using std::vector;
using std::map;
using std::runtime_error;

namespace py = pybind11;

// External function declarations
extern vector<atom_struct> PDBparser(string pdb_filename, string type_filename, bool keep_unknown);
extern void SimpleSolverCL(vector<atom_struct>& pdb, vector<float>& points, int cl_mode);
extern void ChainSelector(vector<vector<string>>& chain_sep, vector<atom_struct>& pdb);
extern void Generic_Solver(vector<atom_struct>& pdb, vector<float>& points, vector<vector<string>> obj1, int mode, int cl_mode);
extern void GeneratePairInteractionData(vector<atom_struct>& pdb);
extern void RelativeSASA(vector<atom_struct>& pdb);
extern void PrintSASAResults(const vector<atom_struct>& pdb, string outname);
extern void PrintSplitAsaAtom(const vector<atom_struct>& pdb, string outname, int mode);
extern void PrintDNA_ProtResults(const vector<atom_struct>& pdb, string outname);
extern void PrintDSASAResults(const vector<atom_struct>& pdb, string outname);
extern void PrintDNA_ProtResultsByAtomMatrix(const vector<atom_struct>& pdb, string basename, int mode);
extern void Print_MatrixInsideAtom(const vector<atom_struct>& pdb, string basename, int mode);


class DrSASA {
private:
    float probe_radius;
    std::unique_ptr<VDWcontainer> vdw_radii;
    std::shared_ptr<std::vector<atom_struct>> atoms;
    bool initialized;
    int cl_mode;
    bool matrix_output;
    bool auto_sort;
    GPUConfig::ComputeMode compute_mode;  

public:
    enum Mode {
        SIMPLE_SASA = 0,         // Simple SASA solver
        CHAIN_DELTA = 1,         // Delta SASA by chain
        RESIDUE_DSASA = 2,       // Residue dSASA mode
        ATOM_DSASA = 3,          // Atom dSASA mode
        CONTACT_SURFACE = 4,     // Atom contact surface mode
        RELATIVE_ASA = 100       // Relative ASA calculation
    };

    enum ComputeBackend {
        CPU = 0,
        OPENCL = 1,
        CUDA = 2
    };

    py::dict calculate_sasa(const string& pdb_file, 
                        bool keep_hetatm = true,
                        const string& output_prefix = "") {
        if (!initialized) {
            throw std::runtime_error("DrSASA not properly initialized");
        }

        try {
            atoms->clear();
            *atoms = PDBparser(pdb_file, "", keep_hetatm);
            if (atoms->empty()) {
                throw std::runtime_error("No atoms read from PDB file");
            }

            if (!vdw_radii) {
                throw std::runtime_error("VDW radii container not initialized");
            }
            vdw_radii->SetRadius(*atoms, probe_radius);

            if (vdw_radii->Points.empty()) {
                throw std::runtime_error("No surface points generated");
            }
            SimpleSolverCL(*atoms, vdw_radii->Points, cl_mode);

    // Simple SASA (Mode 0)
    py::dict calculate_delta_sasa(const string& pdb_file,
                                const vector<vector<string>>& chains,
                                const string& output_prefix) {
        atoms = PDBparser(pdb_file, "", true);
        auto chain_copy = chains;  // Create a non-const copy
        ChainSelector(chain_copy, atoms);  // Pass the copy instead
        vdw_radii.SetRadius(atoms, probe_radius);
        Generic_Solver(atoms, vdw_radii.Points, chains, CHAIN_DELTA, cl_mode);
        GeneratePairInteractionData(atoms);
        
        if (!output_prefix.empty()) {
            write_delta_output(output_prefix);
        }
        
        return create_interface_results_dict();
    }
    // Chain Delta SASA (Mode 1)
    py::dict calculate_residue_sasa(const string& pdb_file,
                                        const vector<string>& chain,
                                        const string& output_prefix) {
        atoms = PDBparser(pdb_file, "", true);
        vector<vector<string>> chain_sep = {chain};
        ChainSelector(chain_sep, atoms);
        vdw_radii.SetRadius(atoms, probe_radius);
        Generic_Solver(atoms, vdw_radii.Points, chain_sep, RESIDUE_DSASA, cl_mode);
        
        if (!output_prefix.empty()) {
            write_output(output_prefix, STANDARD);
        }
        
        return create_results_dict();
    }
    // Residue dSASA (Mode 2)
    py::dict calculate_atom_sasa(const string& pdb_file,
                                        const vector<string>& chain,
                                        const string& output_prefix) {
        atoms = PDBparser(pdb_file, "", true);
        vector<vector<string>> chain_sep = {chain};
        ChainSelector(chain_sep, atoms);
        vdw_radii.SetRadius(atoms, probe_radius);
        Generic_Solver(atoms, vdw_radii.Points, chain_sep, ATOM_DSASA, cl_mode);
        
        if (!output_prefix.empty()) {
            write_output(output_prefix, STANDARD);
        }
        
        return create_results_dict();
    }
    // Atom dSASA (Mode 3)

    // Contact Surface (Mode 4)
    py::dict calculate_contact_surface(const string& pdb_file,
                                            const vector<string>& chain,
                                            const string& output_prefix) {
        atoms = PDBparser(pdb_file, "", true);
        if (!chain.empty()) {
            vector<vector<string>> chain_sep = {chain};
            ChainSelector(chain_sep, atoms);
        }
        vdw_radii.SetRadius(atoms, probe_radius);
        Generic_Solver(atoms, vdw_radii.Points, {chain}, CONTACT_SURFACE, cl_mode);
        
        if (!output_prefix.empty()) {
            write_output(output_prefix, STANDARD);
        }
        
        return create_results_dict();
    }
    // Relative ASA (Mode 100)
    py::dict calculate_relative_asa(const string& pdb_file,
                                        const string& output_prefix) {
        atoms = PDBparser(pdb_file, "", true);
        vdw_radii.SetRadius(atoms, probe_radius);
        SimpleSolverCL(atoms, vdw_radii.Points, cl_mode);
        RelativeSASA(atoms);
        
        if (!output_prefix.empty()) {
            write_output(output_prefix, STANDARD);
        }
        
        return create_relative_results_dict();
    }
    // Simple SASA (Mode 0) py::dict calculate_sasa()
    // Chain Delta SASA (Mode 1) py::dict calculate_delta_sasa();
    // Residue dSASA (Mode 2) py::dict calculate_residue_sasa();
    // Atom dSASA (Mode 3) py::dict calculate_atom_sasa();
    // Contact Surface (Mode 4) py::dict calculate_contact_surface();
    // Relative ASA (Mode 100)  py::dict calculate_relative_asa();

    // Constructor
    DrSASA(const string& vdw_file = "vdw.radii", 
            float probe_radius = 1.4,
            ComputeBackend backend = CPU) 
            : probe_radius(probe_radius),
            initialized(false),
            cl_mode(0),
            matrix_output(true),
            auto_sort(true),
            compute_mode(convert_backend(backend))
        {
            try {
                atoms = std::make_shared<std::vector<atom_struct>>();
                vdw_radii = std::make_unique<VDWcontainer>(vdw_file);
                vdw_radii->GenPoints();
                initialize_compute_device();
                initialized = true;
            } catch (const std::exception& e) {
                throw std::runtime_error("Failed to initialize DrSASA: " + std::string(e.what()));
            }
    }

    private:
        enum OutputMode {
            STANDARD = 0,
            ATMASA_SORT = 1,
            ATMASA_NO_SORT = 2,
            BSA_SORT = 3
        };

        struct ResidueData {
            double sasa;
            vector<const atom_struct*> atoms;
            string chain;
            int number;
            string name;
        };

        struct ContactAnalysis {
            double total_area;
            double polar_area;
            double nonpolar_area;
            vector<pair<int, int>> atom_pairs;
            map<string, double> residue_contributions;
        };

        struct ChainAnalysis {
            string chain_id;
            double total_sasa;
            double buried_sasa;
            double exposed_sasa;
            vector<string> interacting_chains;
            map<string, ContactAnalysis> chain_contacts;
        };

    public:
        // Basic methods
        void set_compute_backend(ComputeBackend backend) {
            compute_mode = convert_backend(backend);
            initialize_compute_device();
        };

        ComputeBackend get_compute_backend() const {
            return convert_mode(compute_mode);
        };

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

        // Configuration methods
        void set_probe_radius(float radius) { probe_radius = radius; }
        float get_probe_radius() const { return probe_radius; }
        void set_matrix_output(bool enable) { matrix_output = enable; }
        bool get_matrix_output() const { return matrix_output; }
        void set_auto_sort(bool enable) { auto_sort = enable; }
        bool get_auto_sort() const { return auto_sort; }
        void set_cl_mode(int mode) { cl_mode = mode; }
        int get_cl_mode() const { return cl_mode; }

public:
    // Advanced chain analysis
    py::dict analyze_chain_interactions(const string& pdb_file,
                                      const vector<string>& chains) {
        atoms = PDBparser(pdb_file, "", true);
        vector<ChainAnalysis> chain_analyses;
        
        // Validate chains
        set<string> available_chains;
        for (const auto& atom : atoms) {
            available_chains.insert(atom.CHAIN);
        }
        
        for (const auto& chain : chains) {
            if (available_chains.find(chain) == available_chains.end()) {
                throw runtime_error("Chain " + chain + " not found in structure");
            }
        }

        // Calculate individual chain SASAs first
        map<string, double> chain_sasas;
        for (const auto& chain : chains) {
            vector<vector<string>> single_chain = {{chain}};
            ChainSelector(single_chain, atoms);
            vdw_radii.SetRadius(atoms, probe_radius);
            SimpleSolverCL(atoms, vdw_radii.Points, cl_mode);
            
            double chain_sasa = 0.0;
            for (const auto& atom : atoms) {
                chain_sasa += atom.SASA;
            }
            chain_sasas[chain] = chain_sasa;
        }

        // Calculate chain-chain interactions
        for (size_t i = 0; i < chains.size(); ++i) {
            ChainAnalysis analysis;
            analysis.chain_id = chains[i];
            analysis.total_sasa = chain_sasas[chains[i]];
            
            for (size_t j = i + 1; j < chains.size(); ++j) {
                vector<vector<string>> chain_pair = {{chains[i]}, {chains[j]}};
                ChainSelector(chain_pair, atoms);
                vdw_radii.SetRadius(atoms, probe_radius);
                Generic_Solver(atoms, vdw_radii.Points, chain_pair, 1, cl_mode);
                GeneratePairInteractionData(atoms);
                
                ContactAnalysis contact;
                analyze_contacts(contact);
                analysis.chain_contacts[chains[j]] = contact;
                analysis.interacting_chains.push_back(chains[j]);
            }
            
            chain_analyses.push_back(analysis);
        }

        return create_chain_analysis_dict(chain_analyses);
    }

    // Advanced residue interaction analysis
    py::dict analyze_residue_interactions(const string& pdb_file,
                                        const string& chain,
                                        double distance_cutoff = 4.5) {
        atoms = PDBparser(pdb_file, "", true);
        
        if (!chain.empty()) {
            vector<vector<string>> chain_sep = {{chain}};
            ChainSelector(chain_sep, atoms);
        }
        
        vdw_radii.SetRadius(atoms, probe_radius);
        Generic_Solver(atoms, vdw_radii.Points, {{chain}}, 2, cl_mode);
        GeneratePairInteractionData(atoms);
        
        map<string, vector<string>> residue_contacts;
        for (const auto& atom : atoms) {
            string res_key = atom.CHAIN + "_" + atom.RESN + "_" + to_string(atom.RESI);
            
            for (const auto& partner_idx : atom.INTERACTION_P) {
                const auto& partner = atoms[partner_idx];
                
                // Calculate distance between atoms
                double dx = atom.COORDS[0] - partner.COORDS[0];
                double dy = atom.COORDS[1] - partner.COORDS[1];
                double dz = atom.COORDS[2] - partner.COORDS[2];
                double distance = sqrt(dx*dx + dy*dy + dz*dz);
                
                // Only consider contacts within distance_cutoff
                if (distance <= distance_cutoff) {
                    string partner_key = partner.CHAIN + "_" + partner.RESN + "_" + 
                                    to_string(partner.RESI);
                    
                    if (res_key != partner_key) {
                        residue_contacts[res_key].push_back(partner_key);
                    }
                }
            }
        }
        
        return create_residue_interaction_dict(residue_contacts);
    }

private:
    void write_output(const string& prefix, OutputMode mode) {
        string asa_file = prefix + ".asa.pdb";
        string atmasa_file = prefix + ".atmasa";
        
        PrintSASAResults(atoms, asa_file);
        PrintSplitAsaAtom(atoms, atmasa_file, static_cast<int>(mode));
        
        if (matrix_output) {
            PrintDNA_ProtResultsByAtomMatrix(atoms, prefix, 0);
            Print_MatrixInsideAtom(atoms, prefix, 0);
        }
    }

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

    void analyze_contacts(ContactAnalysis& contact) {
        contact.total_area = 0.0;
        contact.polar_area = 0.0;
        contact.nonpolar_area = 0.0;
        
        for (const auto& atom : atoms) {
            if (atom.EXT1 > 0) {
                contact.total_area += atom.EXT1;
                
                bool is_polar = atom.ELEMENT == "N" || atom.ELEMENT == "O" || 
                              atom.ELEMENT == "S";
                if (is_polar) {
                    contact.polar_area += atom.EXT1;
                } else {
                    contact.nonpolar_area += atom.EXT1;
                }
                
                string res_key = atom.CHAIN + "_" + atom.RESN + "_" + 
                               to_string(atom.RESI);
                contact.residue_contributions[res_key] += atom.EXT1;
                
                for (const auto& partner_idx : atom.INTERACTION_P) {
                    contact.atom_pairs.emplace_back(atom.ID, atoms[partner_idx].ID);
                }
            }
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

    static GPUConfig::ComputeMode convert_backend(ComputeBackend backend) {
        #ifdef __APPLE__
            return GPUConfig::ComputeMode::CPU;
        #else
            switch(backend) {
                case CPU: return GPUConfig::ComputeMode::CPU;
                case OPENCL: return GPUConfig::ComputeMode::OpenCL;
                case CUDA: return GPUConfig::ComputeMode::CUDA;
                default: return GPUConfig::ComputeMode::CPU;
            }
        #endif
    }

    static ComputeBackend convert_mode(GPUConfig::ComputeMode mode) {
        switch(mode) {
            case GPUConfig::ComputeMode::CPU: return CPU;
            case GPUConfig::ComputeMode::CUDA: return CUDA;
            default: return CPU;
        }
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

    py::dict create_chain_analysis_dict(const vector<ChainAnalysis>& analyses) const {
        py::dict result;
        py::list chain_list;
        
        for (const auto& analysis : analyses) {
            py::dict chain_dict;
            chain_dict["chain_id"] = analysis.chain_id;
            chain_dict["total_sasa"] = analysis.total_sasa;
            chain_dict["buried_sasa"] = analysis.buried_sasa;
            chain_dict["exposed_sasa"] = analysis.exposed_sasa;
            
            py::dict contacts;
            for (const auto& [chain_id, contact] : analysis.chain_contacts) {
                py::dict contact_dict;
                contact_dict["total_area"] = contact.total_area;
                contact_dict["polar_area"] = contact.polar_area;
                contact_dict["nonpolar_area"] = contact.nonpolar_area;
                
                py::dict residue_contribs;
                for (const auto& [res_key, area] : contact.residue_contributions) {
                    residue_contribs[res_key.c_str()] = area;
                }
                contact_dict["residue_contributions"] = residue_contribs;
                
                contacts[chain_id.c_str()] = contact_dict;
            }
            chain_dict["contacts"] = contacts;
            
            chain_list.append(chain_dict);
        }
        
        result["chains"] = chain_list;
        return result;
    }

    py::dict create_residue_interaction_dict(const map<string, vector<string>>& contacts) const {
        py::dict results;
        py::dict residue_contacts;
        
        for (const auto& [res_key, partners] : contacts) {
            py::list partner_list;
            for (const auto& partner : partners) {
                partner_list.append(partner);
            }
            residue_contacts[res_key.c_str()] = partner_list;
        }
        
        results["residue_contacts"] = residue_contacts;
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
};


PYBIND11_MODULE(dr_sasa_py, m) {
    m.doc() = "Python bindings for dr_sasa: High-level interface";
    py::register_exception<std::runtime_error>(m, "DrSASAError");
    // Enum definition
    py::enum_<DrSASA::ComputeBackend>(m, "ComputeBackend")
        .value("CPU", DrSASA::CPU)
        .value("OPENCL", DrSASA::OPENCL)
        .value("CUDA", DrSASA::CUDA)
        .export_values();

    py::enum_<DrSASA::Mode>(m, "Mode")
        .value("SIMPLE_SASA", DrSASA::SIMPLE_SASA)
        .value("CHAIN_DELTA", DrSASA::CHAIN_DELTA)
        .value("RESIDUE_DSASA", DrSASA::RESIDUE_DSASA)
        .value("ATOM_DSASA", DrSASA::ATOM_DSASA)
        .value("CONTACT_SURFACE", DrSASA::CONTACT_SURFACE)
        .value("RELATIVE_ASA", DrSASA::RELATIVE_ASA)
        .export_values();
    // Class definition
    py::class_<DrSASA>(m, "DrSASA")
        // Constructor
        .def(py::init<const string&, float, DrSASA::ComputeBackend>(),
             py::arg("vdw_file") = "vdw.radii",
             py::arg("probe_radius") = 1.4,
             py::arg("backend") = DrSASA::CPU)
        
        // Basic calculation methods
        .def("calculate_sasa", &DrSASA::calculate_sasa,
             py::arg("pdb_file"),
             py::arg("keep_hetatm") = true,
             py::arg("output_prefix") = "")
        .def("calculate_delta_sasa", &DrSASA::calculate_delta_sasa,
             py::arg("pdb_file"),
             py::arg("chains"),
             py::arg("output_prefix") = "")
             
        // Advanced analysis methods
        .def("analyze_chain_interactions", &DrSASA::analyze_chain_interactions,
             py::arg("pdb_file"),
             py::arg("chains"))
        .def("analyze_residue_interactions", &DrSASA::analyze_residue_interactions,
             py::arg("pdb_file"),
             py::arg("chain"),
             py::arg("distance_cutoff") = 4.5)
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

        // Configuration methods
        .def("set_compute_backend", &DrSASA::set_compute_backend)
        .def("get_compute_backend", &DrSASA::get_compute_backend)
        .def("get_compute_device_info", &DrSASA::get_compute_device_info)
        
        // Properties
        .def_property("probe_radius",
            [](const DrSASA& self) { return self.get_probe_radius(); },
            [](DrSASA& self, float radius) {
                if (radius <= 0.0f) throw std::invalid_argument("Probe radius must be positive");
                self.set_probe_radius(radius);
            })
        // Add context manager support
        .def("__enter__", [](DrSASA& self) { return &self; })
        .def("__exit__", [](DrSASA& self, py::object type, py::object value, py::object traceback) {})  
        .def_property("matrix_output",
                     &DrSASA::get_matrix_output,
                     &DrSASA::set_matrix_output)
        .def_property("auto_sort",
                     &DrSASA::get_auto_sort,
                     &DrSASA::set_auto_sort)
        .def_property("cl_mode",
                     &DrSASA::get_cl_mode,
                     &DrSASA::set_cl_mode);
}
