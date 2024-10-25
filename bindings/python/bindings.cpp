// Standard library includes
#include <string>
#include <vector>
#include <cstdint>
#include <map>
#include <sstream>
#include <stdexcept>
#include <set>
#include <cmath>

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

namespace py = pybind11;

using std::string;
using std::vector;
using std::map;
using std::set;
using std::pair;
using std::runtime_error;
using std::to_string;

// External function declarations
extern void SimpleSolverCL(vector<atom_struct>& pdb,  vector<float>& points, int cl_mode);
extern void ChainSelector(vector<vector<string>>& chain_sep,  vector<atom_struct>& pdb);
extern void Generic_Solver(vector<atom_struct>& pdb, vector<float>& points, vector<vector<string>> obj1, int mode, int cl_mode);
extern void RelativeSASA(vector<atom_struct>& pdb);

namespace Constants {
    constexpr double DEFAULT_PROBE_RADIUS = 1.4;
    constexpr double DEFAULT_DISTANCE_CUTOFF = 4.5;
    constexpr const char* DEFAULT_VDW_FILE = "vdw.radii";
    constexpr const char* DEFAULT_OUTPUT_PREFIX = "";
}

class DrSASA {
public:
    // Calculation modes enum
    enum class Mode : int {
        SIMPLE_SASA = 0,         // Simple SASA solver
        CHAIN_DELTA = 1,         // Delta SASA by chain
        RESIDUE_DSASA = 2,       // Residue dSASA mode
        ATOM_DSASA = 3,          // Atom dSASA mode
        CONTACT_SURFACE = 4,     // Atom contact surface mode
        RELATIVE_ASA = 100       // Relative ASA calculation
    };

    // Compute backend enum
    enum class ComputeBackend : int {
        CPU = 0,
        OPENCL = 1,
        CUDA = 2
    };

private:
    // Internal calculation modes
    enum class OutputMode : int {
        STANDARD = 0,
        ATMASA_SORT = 1,
        ATMASA_NO_SORT = 2,
        BSA_SORT = 3
    };

    // Data structures for analysis results
    struct ResidueData {
        double sasa{0.0};
        vector<const atom_struct*> atoms;
        string chain;
        int number{0};
        string name;
    };

    struct ContactAnalysis {
        double total_area{0.0};
        double polar_area{0.0};
        double nonpolar_area{0.0};
        vector<pair<int, int>> atom_pairs;
        map<string, double> residue_contributions;
    };

    struct ChainAnalysis {
        string chain_id;
        double total_sasa{0.0};
        double buried_sasa{0.0};
        double exposed_sasa{0.0};
        vector<string> interacting_chains;
        map<string, ContactAnalysis> chain_contacts;
    };

    // Member variables
    float probe_radius_;
    VDWcontainer vdw_radii_;
    vector<atom_struct> atoms_;
    int cl_mode_{0};
    bool matrix_output_{true};
    bool auto_sort_{true};
    GPUConfig::ComputeMode compute_mode_;
    bool initialized_{false};

public:
    // Constructor and destructor
    explicit DrSASA(const string& vdw_file = "vdw.radii", 
                   float probe_radius = 1.4f,
                   ComputeBackend backend = ComputeBackend::CPU);
    ~DrSASA();

    // Delete copy constructor and assignment operator
    DrSASA(const DrSASA&) = delete;
    DrSASA& operator=(const DrSASA&) = delete;

    // Allow move operations
    DrSASA(DrSASA&&) noexcept = default;
    DrSASA& operator=(DrSASA&&) noexcept = default;

    // Prevent copying and moving after initialization
    void ensure_not_initialized() const {
        if (initialized_) {
            throw std::runtime_error("Cannot modify DrSASA after initialization");
        }
    }
    // Core calculation methods
    py::dict calculate_sasa(const string& pdb_file, 
                           bool keep_hetatm = true,
                           const string& output_prefix = "") {
        if (!initialized_) {
            throw std::runtime_error("DrSASA not properly initialized");
        }
        
        try {
            atoms_ = PDBparser(pdb_file, "", keep_hetatm);
            if (atoms_.empty()) {
                throw std::runtime_error("No atoms loaded from PDB file");
            }
            
            vdw_radii_.SetRadius(atoms_, probe_radius_);
            SimpleSolverCL(atoms_, vdw_radii_.Points, cl_mode_);
            
            return create_results_dict();
        } catch (const std::exception& e) {
            throw std::runtime_error(std::string("SASA calculation failed: ") + e.what());
        }
    }

    py::dict calculate_delta_sasa(const string& pdb_file,
                                 const vector<vector<string>>& chains,
                                 const string& output_prefix = "") {
        if (chains.empty()) {
            throw std::invalid_argument("No chains specified for delta SASA calculation");
        }
        
        try {
            atoms_ = PDBparser(pdb_file, "", true);
            auto chain_copy = chains;
            ChainSelector(chain_copy, atoms_);
            
            if (atoms_.empty()) {
                throw std::runtime_error("No atoms remaining after chain selection");
            }
            
            vdw_radii_.SetRadius(atoms_, probe_radius_);
            Generic_Solver(atoms_, vdw_radii_.Points, chains, 
                         static_cast<int>(Mode::CHAIN_DELTA), cl_mode_);
            GeneratePairInteractionData(atoms_);
            
            return create_interface_results_dict();
        } catch (const std::exception& e) {
            throw std::runtime_error(std::string("Delta SASA calculation failed: ") + e.what());
        }
    }

    // Advanced analysis methods
    py::dict analyze_chain_interactions(const string& pdb_file,
                                      const vector<string>& chains,
                                      double cutoff = 4.5) {
        validate_chains(chains);
        try {
            atoms_ = PDBparser(pdb_file, "", true);
            vector<ChainAnalysis> chain_analyses = perform_chain_analysis(chains, cutoff);
            return create_chain_analysis_dict(chain_analyses);
        } catch (const std::exception& e) {
            throw std::runtime_error(std::string("Chain interaction analysis failed: ") + e.what());
        }
    }

    py::dict analyze_residue_interactions(const string& pdb_file,
                                        const string& chain,
                                        double distance_cutoff = 4.5) {
        try {
            atoms_ = PDBparser(pdb_file, "", true);
            if (!chain.empty()) {
                vector<vector<string>> chain_sep = {{chain}};
                ChainSelector(chain_sep, atoms_);
            }
            
            return perform_residue_analysis(distance_cutoff);
        } catch (const std::exception& e) {
            throw std::runtime_error(std::string("Residue interaction analysis failed: ") + e.what());
        }
    }

    // Configuration methods
    void set_compute_backend(ComputeBackend backend) {
        ensure_not_initialized();
        compute_mode_ = convert_backend(backend);
        initialize_compute_device();
    }

    ComputeBackend get_compute_backend() const {
        return convert_mode(compute_mode_);
    }

    string get_compute_device_info() const {
        switch(compute_mode_) {
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

    // Property getters and setters
    void set_probe_radius(float radius) {
        ensure_not_initialized();
        if (radius <= 0.0f) {
            throw std::invalid_argument("Probe radius must be positive");
        }
        probe_radius_ = radius;
    }

    float get_probe_radius() const { return probe_radius_; }
    
    void set_matrix_output(bool enable) { matrix_output_ = enable; }
    bool get_matrix_output() const { return matrix_output_; }
    
    void set_auto_sort(bool enable) { auto_sort_ = enable; }
    bool get_auto_sort() const { return auto_sort_; }
    
    void set_cl_mode(int mode) {
        ensure_not_initialized();
        cl_mode_ = mode;
    }
    
    int get_cl_mode() const { return cl_mode_; }

private:
    // Initialization and validation methods
    void initialize_compute_device() {
        #ifdef __APPLE__
            compute_mode_ = GPUConfig::ComputeMode::CPU;
            cl_mode_ = 0;
        #else
            switch(compute_mode_) {
                case GPUConfig::ComputeMode::CUDA:
                    if (!check_cuda_availability()) {
                        throw std::runtime_error("CUDA device not available");
                    }
                    cl_mode_ = 1;
                    break;
                case GPUConfig::ComputeMode::OpenCL:
                    if (!check_opencl_availability()) {
                        throw std::runtime_error("OpenCL device not available");
                    }
                    cl_mode_ = 2;
                    break;
                default:
                    cl_mode_ = 0;
                    break;
            }
        #endif
        initialized_ = true;
    }

    bool check_cuda_availability() const {
        // Implementation depends on CUDA runtime
        #ifdef USE_CUDA
            int deviceCount = 0;
            cudaGetDeviceCount(&deviceCount);
            return deviceCount > 0;
        #else
            return false;
        #endif
    }

    bool check_opencl_availability() const {
        // Implementation depends on OpenCL runtime
        #ifdef USE_OPENCL
            // OpenCL availability check implementation
            return true; // Placeholder
        #else
            return false;
        #endif
    }

    void validate_chains(const vector<string>& chains) const {
        if (chains.empty()) {
            throw std::invalid_argument("Empty chain list provided");
        }

        set<string> available_chains;
        for (const auto& atom : atoms_) {
            available_chains.insert(atom.CHAIN);
        }

        for (const auto& chain : chains) {
            if (available_chains.find(chain) == available_chains.end()) {
                throw std::invalid_argument("Chain " + chain + " not found in structure");
            }
        }
    }

    // Analysis helper methods
    vector<ChainAnalysis> perform_chain_analysis(const vector<string>& chains, 
                                               double cutoff) {
        vector<ChainAnalysis> analyses;
        map<string, double> chain_sasas = calculate_individual_chain_sasas(chains);

        for (size_t i = 0; i < chains.size(); ++i) {
            ChainAnalysis analysis;
            analysis.chain_id = chains[i];
            analysis.total_sasa = chain_sasas[chains[i]];

            for (size_t j = i + 1; j < chains.size(); ++j) {
                auto contact = analyze_chain_pair(chains[i], chains[j], cutoff);
                if (contact.total_area > 0.0) {
                    analysis.chain_contacts[chains[j]] = contact;
                    analysis.interacting_chains.push_back(chains[j]);
                }
            }

            analyses.push_back(std::move(analysis));
        }

        return analyses;
    }

    map<string, double> calculate_individual_chain_sasas(
        const vector<string>& chains) {
        map<string, double> chain_sasas;
        
        for (const auto& chain : chains) {
            vector<vector<string>> single_chain = {{chain}};
            vector<atom_struct> chain_atoms = atoms_; // Make a copy
            
            ChainSelector(single_chain, chain_atoms);
            vdw_radii_.SetRadius(chain_atoms, probe_radius_);
            SimpleSolverCL(chain_atoms, vdw_radii_.Points, cl_mode_);
            
            double chain_sasa = 0.0;
            for (const auto& atom : chain_atoms) {
                chain_sasa += atom.SASA;
            }
            chain_sasas[chain] = chain_sasa;
        }
        
        return chain_sasas;
    }

    ContactAnalysis analyze_chain_pair(const string& chain1, 
                                     const string& chain2,
                                     double cutoff) {
        vector<vector<string>> chain_pair = {{chain1}, {chain2}};
        vector<atom_struct> pair_atoms = atoms_; // Make a copy
        
        ChainSelector(chain_pair, pair_atoms);
        vdw_radii_.SetRadius(pair_atoms, probe_radius_);
        Generic_Solver(pair_atoms, vdw_radii_.Points, chain_pair, 
                      static_cast<int>(Mode::CHAIN_DELTA), cl_mode_);
        GeneratePairInteractionData(pair_atoms);
        
        ContactAnalysis contact;
        analyze_contacts(pair_atoms, contact, cutoff);
        return contact;
    }

    void process_contact_interactions(const atom_struct& atom,
                                    const vector<atom_struct>& atoms,
                                    ContactAnalysis& contact,
                                    double cutoff) const {
        for (const auto& partner_idx : atom.INTERACTION_P) {
            const auto& partner = atoms[partner_idx];
            
            float dx = atom.COORDS[0] - partner.COORDS[0];
            float dy = atom.COORDS[1] - partner.COORDS[1];
            float dz = atom.COORDS[2] - partner.COORDS[2];
            float distance = std::sqrt(dx*dx + dy*dy + dz*dz);

            if (distance <= cutoff) {
                contact.atom_pairs.emplace_back(atom.ID, partner.ID);
            }
        }
    }


    void analyze_contacts(const vector<atom_struct>& atoms,
                         ContactAnalysis& contact,
                         double cutoff) const {
        for (const auto& atom : atoms) {
            if (atom.EXT1 > 0) {
                process_contact_interactions(atom, atoms, contact, cutoff);
            }
        }
    }

    static void update_contact_areas(const atom_struct& atom,
                                   ContactAnalysis& contact) {
        contact.total_area += atom.EXT1;
        
        bool is_polar = atom.ELEMENT == "N" || 
                       atom.ELEMENT == "O" || 
                       atom.ELEMENT == "S";
        
        if (is_polar) {
            contact.polar_area += atom.EXT1;
        } else {
            contact.nonpolar_area += atom.EXT1;
        }
        
        string res_key = atom.CHAIN + "_" + 
                        atom.RESN + "_" + 
                        to_string(atom.RESI);
        contact.residue_contributions[res_key] += atom.EXT1;
    }

    py::dict perform_residue_analysis(double distance_cutoff) {
        try {
            vdw_radii_.SetRadius(atoms_, probe_radius_);
            Generic_Solver(atoms_, vdw_radii_.Points, {}, 
                         static_cast<int>(Mode::RESIDUE_DSASA), cl_mode_);
            GeneratePairInteractionData(atoms_);

            map<string, vector<string>> residue_contacts;
            process_atom_interactions(atoms_, residue_contacts, distance_cutoff);

            return create_residue_interaction_dict(residue_contacts);
        } catch (const std::exception& e) {
            throw std::runtime_error("Residue analysis failed: " + std::string(e.what()));
        }
    }

    void process_atom_interactions(const vector<atom_struct>& atoms,
                                 map<string, vector<string>>& residue_contacts,
                                 double distance_cutoff) const {
        for (const auto& atom : atoms) {
            string res_key = atom.CHAIN + "_" + atom.RESN + "_" + 
                           std::to_string(atom.RESI);

            for (const auto& partner_idx : atom.INTERACTION_P) {
                const auto& partner = atoms[partner_idx];
                
                // Calculate distance between atoms
                float dx = atom.COORDS[0] - partner.COORDS[0];
                float dy = atom.COORDS[1] - partner.COORDS[1];
                float dz = atom.COORDS[2] - partner.COORDS[2];
                float distance = std::sqrt(dx*dx + dy*dy + dz*dz);

                if (distance <= distance_cutoff) {
                    string partner_key = partner.CHAIN + "_" + 
                                       partner.RESN + "_" + 
                                       std::to_string(partner.RESI);

                    if (res_key != partner_key) {
                        residue_contacts[res_key].push_back(partner_key);
                    }
                }
            }
        }
    }

    // Result dictionary creation methods
    py::dict create_results_dict() const {
        py::dict results;
        try {
            // Basic atom information
            vector<double> atom_sasa;
            vector<string> atom_names;
            vector<string> residue_names;
            vector<int> residue_numbers;
            vector<string> chain_ids;
            vector<string> elements;
            vector<float> coordinates;
            
            atom_sasa.reserve(atoms_.size());
            atom_names.reserve(atoms_.size());
            residue_names.reserve(atoms_.size());
            residue_numbers.reserve(atoms_.size());
            chain_ids.reserve(atoms_.size());
            elements.reserve(atoms_.size());
            coordinates.reserve(atoms_.size() * 3);
            
            double total_sasa = 0.0;
            
            // Populate vectors
            for (const auto& atom : atoms_) {
                atom_sasa.push_back(atom.SASA);
                atom_names.push_back(atom.NAME);
                residue_names.push_back(atom.RESN);
                residue_numbers.push_back(atom.RESI);
                chain_ids.push_back(atom.CHAIN);
                elements.push_back(atom.ELEMENT);
                
                // Coordinates as flat array
                coordinates.push_back(atom.COORDS[0]);
                coordinates.push_back(atom.COORDS[1]);
                coordinates.push_back(atom.COORDS[2]);
                
                total_sasa += atom.SASA;
            }
            
            // Create numpy arrays for better Python integration
            results["atom_sasa"] = py::array_t<double>(atom_sasa.size(), atom_sasa.data());
            results["atom_names"] = atom_names;
            results["residue_names"] = residue_names;
            results["residue_numbers"] = py::array_t<int>(residue_numbers.size(), residue_numbers.data());
            results["chain_ids"] = chain_ids;
            results["elements"] = elements;
            results["coordinates"] = create_coordinates_array();
            results["total_sasa"] = total_sasa;
            
            return results;
        } catch (const std::exception& e) {
            throw std::runtime_error("Failed to create results dictionary: " + std::string(e.what()));
        }
    }

    py::dict create_relative_results_dict() const {
        try {
            py::dict results = create_results_dict();  // Get basic SASA results
            
            vector<double> rel_sasa;
            rel_sasa.reserve(atoms_.size());
            
            for (const auto& atom : atoms_) {
                double relative = (atom.AREA > 0.0) ? 
                    (atom.SASA / atom.AREA * 100.0) : 0.0;  // as percentage
                rel_sasa.push_back(relative);
            }
            
            results["relative_sasa"] = py::array_t<double>(rel_sasa.size(), rel_sasa.data());
            return results;
        } catch (const std::exception& e) {
            throw std::runtime_error("Failed to create relative results dictionary: " + 
                                   std::string(e.what()));
        }
    }

    py::array_t<float> create_coordinates_array() const {
        std::vector<float> coordinates;
        coordinates.reserve(atoms_.size() * 3);

        for (const auto& atom : atoms_) {
            coordinates.push_back(atom.COORDS[0]);
            coordinates.push_back(atom.COORDS[1]);
            coordinates.push_back(atom.COORDS[2]);
        }

        std::vector<ssize_t> shape = {static_cast<ssize_t>(atoms_.size()), 3};
        std::vector<ssize_t> strides = {3 * sizeof(float), sizeof(float)};

        return py::array_t<float>(
            shape,
            strides,
            coordinates.data()
        );
    }

    py::dict create_interface_results_dict() const {
        try {
            py::dict results;
            py::list interface_atoms;
            double total_interface_area = 0.0;

            for (const auto& atom : atoms_) {
                if (atom.EXT1 > 0) {  // EXT1 stores dSASA
                    py::dict atom_info;
                    atom_info["id"] = atom.ID;
                    atom_info["name"] = atom.NAME;
                    atom_info["resname"] = atom.RESN;
                    atom_info["chain"] = atom.CHAIN;
                    atom_info["resid"] = atom.RESI;
                    atom_info["delta_sasa"] = atom.EXT1;
                    atom_info["coordinates"] = py::make_tuple(
                        atom.COORDS[0], atom.COORDS[1], atom.COORDS[2]);
                    
                    interface_atoms.append(atom_info);
                    total_interface_area += atom.EXT1;
                }
            }

            results["interface_atoms"] = interface_atoms;
            results["total_interface_area"] = total_interface_area;
            
            // Add interface statistics
            add_interface_statistics(results);
            
            return results;
        } catch (const std::exception& e) {
            throw std::runtime_error("Failed to create interface results dictionary: " + 
                                   std::string(e.what()));
        }
    }

    void add_interface_statistics(py::dict& results) const {
        double polar_area = 0.0;
        double nonpolar_area = 0.0;
        map<string, double> residue_contributions;
        
        for (const auto& atom : atoms_) {
            if (atom.EXT1 > 0) {
                bool is_polar = atom.ELEMENT == "N" || 
                              atom.ELEMENT == "O" || 
                              atom.ELEMENT == "S";
                
                if (is_polar) {
                    polar_area += atom.EXT1;
                } else {
                    nonpolar_area += atom.EXT1;
                }
                
                string res_key = atom.CHAIN + "_" + 
                               atom.RESN + "_" + 
                               std::to_string(atom.RESI);
                residue_contributions[res_key] += atom.EXT1;
            }
        }
        
        results["polar_area"] = polar_area;
        results["nonpolar_area"] = nonpolar_area;
        results["residue_contributions"] = residue_contributions;
    }

    py::dict create_chain_analysis_dict(const vector<ChainAnalysis>& analyses) const {
        try {
            py::dict result;
            py::list chain_list;
            
            for (const auto& analysis : analyses) {
                py::dict chain_dict;
                chain_dict["chain_id"] = analysis.chain_id;
                chain_dict["total_sasa"] = analysis.total_sasa;
                chain_dict["buried_sasa"] = analysis.buried_sasa;
                chain_dict["exposed_sasa"] = analysis.exposed_sasa;
                
                // Create contacts dictionary
                py::dict contacts;
                for (const auto& [chain_id, contact] : analysis.chain_contacts) {
                    contacts[chain_id.c_str()] = create_contact_dict(contact);
                }
                
                chain_dict["contacts"] = contacts;
                chain_dict["interacting_chains"] = analysis.interacting_chains;
                
                chain_list.append(chain_dict);
            }
            
            result["chains"] = chain_list;
            return result;
        } catch (const std::exception& e) {
            throw std::runtime_error("Failed to create chain analysis dictionary: " + 
                                   std::string(e.what()));
        }
    }

    py::dict create_contact_dict(const ContactAnalysis& contact) const {
        py::dict contact_dict;
        contact_dict["total_area"] = contact.total_area;
        contact_dict["polar_area"] = contact.polar_area;
        contact_dict["nonpolar_area"] = contact.nonpolar_area;
        
        // Convert residue contributions to Python dictionary
        py::dict residue_contribs;
        for (const auto& [res_key, area] : contact.residue_contributions) {
            residue_contribs[res_key.c_str()] = area;
        }
        contact_dict["residue_contributions"] = residue_contribs;
        
        // Add atom pairs if available
        if (!contact.atom_pairs.empty()) {
            py::list pairs;
            for (const auto& [atom1, atom2] : contact.atom_pairs) {
                pairs.append(py::make_tuple(atom1, atom2));
            }
            contact_dict["atom_pairs"] = pairs;
        }
        
        return contact_dict;
    }

    py::dict create_residue_interaction_dict(const map<string, vector<string>>& contacts) const {
        try {
            py::dict results;
            py::dict residue_contacts;
            
            for (const auto& [res_key, partners] : contacts) {
                // Convert vector of partners to Python list
                py::list partner_list;
                for (const auto& partner : partners) {
                    partner_list.append(partner);
                }
                residue_contacts[res_key.c_str()] = partner_list;
            }
            
            results["residue_contacts"] = residue_contacts;
            
            // Add additional residue-level statistics
            add_residue_statistics(results);
            
            return results;
        } catch (const std::exception& e) {
            throw std::runtime_error("Failed to create residue interaction dictionary: " + 
                                   std::string(e.what()));
        }
    }

    void add_residue_statistics(py::dict& results) const {
        try {
            map<string, ResidueData> residue_stats;
            
            // Collect residue-level data
            for (const auto& atom : atoms_) {
                string res_key = atom.CHAIN + "_" + 
                               atom.RESN + "_" + 
                               std::to_string(atom.RESI);
                
                auto& res_data = residue_stats[res_key];
                res_data.sasa += atom.SASA;
                res_data.atoms.push_back(&atom);
                res_data.chain = atom.CHAIN;
                res_data.number = atom.RESI;
                res_data.name = atom.RESN;
            }
            
            // Convert to Python dictionary
            py::dict residue_info;
            for (const auto& [res_key, data] : residue_stats) {
                py::dict res_dict;
                res_dict["sasa"] = data.sasa;
                res_dict["chain"] = data.chain;
                res_dict["number"] = data.number;
                res_dict["name"] = data.name;
                res_dict["n_atoms"] = static_cast<int>(data.atoms.size());
                
                residue_info[res_key.c_str()] = res_dict;
            }
            
            results["residue_info"] = residue_info;
            
        } catch (const std::exception& e) {
            throw std::runtime_error("Failed to add residue statistics: " + 
                                   std::string(e.what()));
        }
    }
    
    static GPUConfig::ComputeMode convert_backend(ComputeBackend backend) {
        #ifdef __APPLE__
            return GPUConfig::ComputeMode::CPU;
        #else
            switch(backend) {
                case ComputeBackend::CPU: 
                    return GPUConfig::ComputeMode::CPU;
                case ComputeBackend::OPENCL: 
                    return GPUConfig::ComputeMode::OpenCL;
                case ComputeBackend::CUDA: 
                    return GPUConfig::ComputeMode::CUDA;
                default: 
                    return GPUConfig::ComputeMode::CPU;
            }
        #endif
    }

    static ComputeBackend convert_mode(GPUConfig::ComputeMode mode) {
        switch(mode) {
            case GPUConfig::ComputeMode::CPU: 
                return ComputeBackend::CPU;
            case GPUConfig::ComputeMode::CUDA: 
                return ComputeBackend::CUDA;
            case GPUConfig::ComputeMode::OpenCL: 
                return ComputeBackend::OPENCL;
            default: 
                return ComputeBackend::CPU;
        }
    }
};

PYBIND11_MODULE(dr_sasa_py, m) {
    m.doc() = R"pbdoc(
        Python bindings for dr_sasa
        ---------------------------

        A high-performance library for calculating Solvent Accessible Surface Area (SASA)
        and analyzing molecular interfaces.

        Key features:
        * SASA calculation
        * Delta SASA analysis
        * Chain interaction analysis
        * Residue interaction analysis
        * Multiple computation backends (CPU, OpenCL, CUDA)
    )pbdoc";

    // Enums
    py::enum_<DrSASA::ComputeBackend>(m, "ComputeBackend", py::arithmetic())
        .value("CPU", DrSASA::ComputeBackend::CPU, "CPU computation using OpenMP")
        .value("OPENCL", DrSASA::ComputeBackend::OPENCL, "GPU computation using OpenCL")
        .value("CUDA", DrSASA::ComputeBackend::CUDA, "GPU computation using CUDA")
        .export_values();

    py::enum_<DrSASA::Mode>(m, "Mode", py::arithmetic())
        .value("SIMPLE_SASA", DrSASA::Mode::SIMPLE_SASA, "Simple SASA calculation")
        .value("CHAIN_DELTA", DrSASA::Mode::CHAIN_DELTA, "Chain-based delta SASA")
        .value("RESIDUE_DSASA", DrSASA::Mode::RESIDUE_DSASA, "Residue-based dSASA")
        .value("ATOM_DSASA", DrSASA::Mode::ATOM_DSASA, "Atom-based dSASA")
        .value("CONTACT_SURFACE", DrSASA::Mode::CONTACT_SURFACE, "Contact surface analysis")
        .value("RELATIVE_ASA", DrSASA::Mode::RELATIVE_ASA, "Relative ASA calculation")
        .export_values();

    // Main class
    py::class_<DrSASA>(m, "DrSASA")
        .def(py::init<const string&, float, DrSASA::ComputeBackend>(),
             py::arg("vdw_file") = "vdw.radii",
             py::arg("probe_radius") = 1.4f,
             py::arg("backend") = DrSASA::ComputeBackend::CPU,
             R"pbdoc(
                Initialize DrSASA calculator.

                Parameters:
                    vdw_file (str): Path to van der Waals radii file
                    probe_radius (float): Probe radius for SASA calculation
                    backend (ComputeBackend): Computation backend to use
             )pbdoc")
        
        // Core calculation methods
        .def("calculate_sasa", &DrSASA::calculate_sasa,
             py::arg("pdb_file"),
             py::arg("keep_hetatm") = true,
             py::arg("output_prefix") = "",
             "Calculate SASA for all atoms in the structure")
        
        .def("calculate_delta_sasa", &DrSASA::calculate_delta_sasa,
             py::arg("pdb_file"),
             py::arg("chains"),
             py::arg("output_prefix") = "",
             "Calculate delta SASA between chains")
        
        // Advanced analysis methods
        .def("analyze_chain_interactions", &DrSASA::analyze_chain_interactions,
            py::arg("pdb_file"),
            py::arg("chains"),
            py::arg("cutoff") = 4.5,
            "Analyze interactions between specified chains")
        
        .def("analyze_residue_interactions", &DrSASA::analyze_residue_interactions,
             py::arg("pdb_file"),
             py::arg("chain"),
             py::arg("distance_cutoff") = 4.5,
             "Analyze residue-level interactions")
        
        // Configuration methods
        .def("set_compute_backend", &DrSASA::set_compute_backend)
        .def("get_compute_backend", &DrSASA::get_compute_backend)
        .def("get_compute_device_info", &DrSASA::get_compute_device_info)

        // Properties
        .def_property("probe_radius",  &DrSASA::get_probe_radius,  &DrSASA::set_probe_radius)
        .def_property("matrix_output", &DrSASA::get_matrix_output, &DrSASA::set_matrix_output)
        .def_property("auto_sort", &DrSASA::get_auto_sort, &DrSASA::set_auto_sort)
        .def_property("cl_mode", &DrSASA::get_cl_mode,  &DrSASA::set_cl_mode);

    // Version information
    m.attr("__version__") = "0.1.0";
}