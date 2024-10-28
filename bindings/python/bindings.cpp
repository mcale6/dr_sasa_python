// base_drsasa.hpp
#pragma once
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <stdexcept>
#include <cstdint>
#include <algorithm>
#include <cmath>
#include <filesystem>

#include "stdafx.h"
#include "atom_struct.h"
#include "PDBparser2.h"
#include "SetRadius.h"
#include "SurfaceSolverCL.h"

namespace py = pybind11;

// Result structure for unified return types
struct DRSASAResult {
    py::array_t<double> sasa;
    py::array_t<double> rel_sasa;
    py::array_t<float> coordinates;
    py::dict metadata;
    py::dict statistics;
    
    static DRSASAResult from_atom_struct(const vector<atom_struct>& atoms);
};

class BaseDRSASA {
protected:
    double probe = 1.4;
    bool keepunknown = true;
    int cl_mode = 0;
    bool reorder = false;
    string vdwfile;
    vector<atom_struct> pdb;
    vector<vector<string>> chain_sep;
    VDWcontainer rad;

    virtual void initialize_vdw() {
        rad = VDWcontainer(vdwfile);
        rad.GenPoints();
        rad.SetRadius(pdb, probe);
    }

    virtual void process_chains() {
        if(chain_sep.size() >= 1) {
            ChainSelector(chain_sep, pdb);
        }
        if(reorder) {
            pdb = ReorderPDB(pdb);
        }
    }

    virtual void load_structure(const string& input_file, const string& types = "") {
        pdb = PDBparser(input_file, types, keepunknown);
    }

    // Pure virtual functions
    virtual void solve() = 0;
    virtual DRSASAResult process_results() = 0;

public:
    virtual ~BaseDRSASA() = default;
    
    void set_probe_radius(double radius) { probe = radius; }
    void set_vdw_file(const string& file) { vdwfile = file; }
    void set_chains(const vector<vector<string>>& chains) { chain_sep = chains; }
    void set_cl_mode(int mode) { cl_mode = mode; }
    void enable_reorder(bool enable) { reorder = enable; }
    
    virtual DRSASAResult run(const string& input_file) {
        load_structure(input_file);
        initialize_vdw();
        process_chains();
        solve();
        return process_results();
    }
};

// simple_sasa.hpp
#pragma once
#include "base_drsasa.hpp"

class SimpleSASA : public BaseDRSASA {
protected:
    void solve() override {
        if (chain_sep.size() > 1) {
            throw std::runtime_error("Simple SASA requires single object");
        }
        SolveInteractions(pdb, 0);
        SimpleSolverCL(pdb, rad.Points, cl_mode);
    }

    DRSASAResult process_results() override {
        DRSASAResult result;
        
        // Prepare arrays
        const size_t n_atoms = pdb.size();
        auto sasa = py::array_t<double>(n_atoms);
        auto metadata = py::dict();
        
        // Process atom data
        double total_sasa = 0.0;
        py::buffer_info buf = sasa.request();
        double* ptr = static_cast<double*>(buf.ptr);
        
        for(size_t i = 0; i < n_atoms; ++i) {
            ptr[i] = pdb[i].SASA;
            total_sasa += pdb[i].SASA;
            
            // Collect metadata
            metadata[pdb[i].CHAIN + "_" + pdb[i].RESN + std::to_string(pdb[i].RESI)] = pdb[i].SASA;
        }
        
        result.sasa = sasa;
        result.metadata = metadata;
        result.statistics = py::dict("total_sasa"_a=total_sasa);
        
        return result;
    }

public:
    SimpleSASA() = default;
};

class RelativeSASA : public SimpleSASA {
protected:
    void solve() override {
        // First do normal SASA calculation
        SimpleSASA::solve();
        
        // Then calculate relative SASA
        RelativeSASA(pdb);
    }

    DRSASAResult process_results() override {
        DRSASAResult result = SimpleSASA::process_results();
        
        // Add relative SASA data
        const size_t n_atoms = pdb.size();
        auto rel_sasa = py::array_t<double>(n_atoms);
        py::buffer_info buf = rel_sasa.request();
        double* ptr = static_cast<double*>(buf.ptr);
        
        for(size_t i = 0; i < n_atoms; ++i) {
            ptr[i] = pdb[i].REL_SASA;
        }
        
        result.rel_sasa = rel_sasa;
        
        // Add statistics about relative SASA
        double avg_rel_sasa = 0.0;
        for(const auto& atom : pdb) {
            avg_rel_sasa += atom.REL_SASA;
        }
        avg_rel_sasa /= n_atoms;
        
        result.statistics["average_relative_sasa"] = avg_rel_sasa;
        
        return result;
    }

public:
    RelativeSASA() = default;
};

// Utility functions for mode 106 (Atom Type SASA)
class AtomTypeSASA : public SimpleSASA {
protected:
    string types;

    void load_structure(const string& input_file, const string& types_def) override {
        types = types_def;
        pdb = PDBparser(input_file, types, keepunknown);
    }

    DRSASAResult process_results() override {
        DRSASAResult result = SimpleSASA::process_results();
        
        // Add atom type specific processing
        map<string, double> sasa_by_type;
        for(const auto& atom : pdb) {
            sasa_by_type[atom.TYPE] += atom.SASA;
        }
        
        result.statistics["sasa_by_type"] = sasa_by_type;
        return result;
    }

public:
    AtomTypeSASA() = default;
    
    DRSASAResult run(const string& input_file, const string& types_def) {
        return SimpleSASA::run(input_file);
    }
};


class GenericDSASA : public BaseDRSASA {
protected:
    int Imode;
    bool mtrx = true;

    void determine_Imode() {
        // Initial mode determination based on chain separators
        if (chain_sep.size() <= 1) {
            Imode = 4;  // Automatic
        } else {
            Imode = 1;  // Manual
        }

        // Check for protein chains
        set<string> proteinChains;
        if (chain_sep.size() <= 1) {
            for (const auto& atom: pdb) {
                if (atom.MOL_TYPE == "PROTEIN") {
                    proteinChains.insert(atom.CHAIN);
                } else {
                    proteinChains.clear();
                    break;
                }
            }
        }
        
        // Special case for exactly two protein chains
        if (proteinChains.size() == 2) {
            Imode = 5;
        }
    }

    void solve() override {
        determine_Imode();
        Generic_Solver(pdb, rad.Points, chain_sep, Imode, cl_mode);
        GeneratePairInteractionData(pdb);
    }

    DRSASAResult process_results() override {
        DRSASAResult result;
        
        const size_t n_atoms = pdb.size();
        auto dsasa = py::array_t<double>(n_atoms);
        auto interaction_matrix = py::array_t<double>({n_atoms, n_atoms});
        
        // Fill dSASA values
        {
            py::buffer_info buf = dsasa.request();
            double* ptr = static_cast<double*>(buf.ptr);
            for(size_t i = 0; i < n_atoms; ++i) {
                ptr[i] = pdb[i].dSASA;
            }
        }
        
        // Process interaction data
        py::dict interactions;
        double total_interface = 0.0;
        map<string, double> interface_by_chain;
        map<pair<string, string>, double> chain_interactions;

        for(size_t i = 0; i < n_atoms; ++i) {
            const auto& atom = pdb[i];
            if(atom.dSASA > 0) {
                interface_by_chain[atom.CHAIN] += atom.dSASA;
                total_interface += atom.dSASA;
                
                // Process pair interactions if available
                for(const auto& interaction : atom.interactions) {
                    auto chain_pair = make_pair(atom.CHAIN, pdb[interaction.first].CHAIN);
                    chain_interactions[chain_pair] += interaction.second;
                }
            }
        }

        // Fill interaction matrix if requested
        if(mtrx) {
            py::buffer_info buf = interaction_matrix.request();
            double* ptr = static_cast<double*>(buf.ptr);
            for(size_t i = 0; i < n_atoms; ++i) {
                for(size_t j = 0; j < n_atoms; ++j) {
                    // Find interaction value between atoms i and j
                    double int_value = 0.0;
                    for(const auto& interaction : pdb[i].interactions) {
                        if(interaction.first == j) {
                            int_value = interaction.second;
                            break;
                        }
                    }
                    ptr[i * n_atoms + j] = int_value;
                }
            }
        }

        // Pack results
        result.sasa = dsasa;
        result.metadata = py::dict(
            "interface_by_chain"_a=interface_by_chain,
            "chain_interactions"_a=chain_interactions
        );
        
        result.statistics = py::dict(
            "total_interface"_a=total_interface,
            "interaction_matrix"_a=interaction_matrix,
            "mode"_a=Imode
        );

        return result;
    }

public:
    GenericDSASA() = default;
    
    void set_matrix_output(bool enable) { mtrx = enable; }
    
    // Override run to include more detailed error checking
    DRSASAResult run(const string& input_file) override {
        if(chain_sep.empty()) {
            throw std::runtime_error("No chain separators defined for Generic dSASA calculation");
        }
        return BaseDRSASA::run(input_file);
    }

    // Helper method to get interaction matrix directly
    py::array_t<double> get_interaction_matrix() {
        if(!mtrx) {
            throw std::runtime_error("Matrix output is disabled");
        }
        
        const size_t n_atoms = pdb.size();
        auto matrix = py::array_t<double>({n_atoms, n_atoms});
        
        py::buffer_info buf = matrix.request();
        double* ptr = static_cast<double*>(buf.ptr);
        
        for(size_t i = 0; i < n_atoms; ++i) {
            for(size_t j = 0; j < n_atoms; ++j) {
                // Find interaction value
                double int_value = 0.0;
                for(const auto& interaction : pdb[i].interactions) {
                    if(interaction.first == j) {
                        int_value = interaction.second;
                        break;
                    }
                }
                ptr[i * n_atoms + j] = int_value;
            }
        }
        
        return matrix;
    }
};

class InternalDSASA : public BaseDRSASA {
protected:
    bool mtrx = true;
    int analysis_mode;  // 2 for residue, 3 for atom
    
    void solve() override {
        if (chain_sep.size() > 1) {
            throw std::runtime_error("Internal dSASA requires single object");
        }
        
        Generic_Solver(pdb, rad.Points, chain_sep, analysis_mode, cl_mode);
        GeneratePairInteractionData(pdb);
    }

    DRSASAResult process_results() override {
        DRSASAResult result;
        const size_t n_atoms = pdb.size();
        
        // Create basic arrays
        auto dsasa = py::array_t<double>(n_atoms);
        py::dict interactions;
        
        // Process based on mode
        if (analysis_mode == 2) {  // Residue mode
            process_residue_results(result);
        } else {  // Atom mode
            process_atom_results(result);
        }
        
        return result;
    }

    void process_residue_results(DRSASAResult& result) {
        // Map to store residue-level data
        map<string, double> residue_sasa;
        map<string, vector<size_t>> residue_atoms;
        map<pair<string, string>, double> residue_interactions;
        
        // First pass: collect residue information
        for(size_t i = 0; i < pdb.size(); ++i) {
            const auto& atom = pdb[i];
            string res_id = atom.CHAIN + "_" + atom.RESN + "_" + to_string(atom.RESI);
            residue_sasa[res_id] += atom.dSASA;
            residue_atoms[res_id].push_back(i);
        }
        
        // Second pass: process interactions
        for(const auto& [res_id, atoms] : residue_atoms) {
            for(size_t atom_idx : atoms) {
                const auto& atom = pdb[atom_idx];
                for(const auto& interaction : atom.interactions) {
                    const auto& partner = pdb[interaction.first];
                    string partner_id = partner.CHAIN + "_" + partner.RESN + "_" + to_string(partner.RESI);
                    if(res_id < partner_id) {  // Avoid double counting
                        residue_interactions[{res_id, partner_id}] += interaction.second;
                    }
                }
            }
        }
        
        // Create interaction matrix if requested
        if(mtrx) {
            const size_t n_residues = residue_sasa.size();
            auto matrix = py::array_t<double>({n_residues, n_residues});
            py::buffer_info buf = matrix.request();
            double* ptr = static_cast<double*>(buf.ptr);
            
            // Create residue index mapping
            map<string, size_t> residue_indices;
            size_t idx = 0;
            for(const auto& [res_id, _] : residue_sasa) {
                residue_indices[res_id] = idx++;
            }
            
            // Fill matrix
            for(const auto& [res_pair, interaction] : residue_interactions) {
                const auto& [res1, res2] = res_pair;
                size_t i = residue_indices[res1];
                size_t j = residue_indices[res2];
                ptr[i * n_residues + j] = interaction;
                ptr[j * n_residues + i] = interaction;  // Symmetric
            }
            
            result.statistics["interaction_matrix"] = matrix;
        }
        
        // Pack results
        result.sasa = create_atom_array(pdb, &atom_struct::dSASA);
        result.metadata = py::dict(
            "residue_sasa"_a=residue_sasa,
            "residue_interactions"_a=residue_interactions
        );
    }

    void process_atom_results(DRSASAResult& result) {
        const size_t n_atoms = pdb.size();
        
        // Create atom-level arrays
        result.sasa = create_atom_array(pdb, &atom_struct::dSASA);
        
        // Process atom interactions
        if(mtrx) {
            auto matrix = py::array_t<double>({n_atoms, n_atoms});
            py::buffer_info buf = matrix.request();
            double* ptr = static_cast<double*>(buf.ptr);
            
            for(size_t i = 0; i < n_atoms; ++i) {
                for(const auto& interaction : pdb[i].interactions) {
                    size_t j = interaction.first;
                    double value = interaction.second;
                    ptr[i * n_atoms + j] = value;
                    ptr[j * n_atoms + i] = value;  // Symmetric
                }
            }
            
            result.statistics["interaction_matrix"] = matrix;
        }
        
        // Collect atom interaction statistics
        map<string, double> interaction_by_type;
        double total_interaction = 0.0;
        
        for(const auto& atom : pdb) {
            string type_key = atom.ELEMENT + "_" + atom.TYPE;
            interaction_by_type[type_key] += atom.dSASA;
            total_interaction += atom.dSASA;
        }
        
        result.metadata = py::dict(
            "interaction_by_type"_a=interaction_by_type,
            "total_interaction"_a=total_interaction
        );
    }

    template<typename T>
    py::array_t<double> create_atom_array(const vector<atom_struct>& atoms, T atom_struct::*member) {
        auto array = py::array_t<double>(atoms.size());
        py::buffer_info buf = array.request();
        double* ptr = static_cast<double*>(buf.ptr);
        
        for(size_t i = 0; i < atoms.size(); ++i) {
            ptr[i] = atoms[i].*member;
        }
        
        return array;
    }

public:
    InternalDSASA(bool residue_mode = true) {
        analysis_mode = residue_mode ? 2 : 3;
    }
    
    void set_matrix_output(bool enable) { mtrx = enable; }
    bool is_residue_mode() const { return analysis_mode == 2; }
    
    static InternalDSASA create_residue_mode() {
        return InternalDSASA(true);
    }
    
    static InternalDSASA create_atom_mode() {
        return InternalDSASA(false);
    }
};

// decoupled_dsasa.hpp
#pragma once
#include "base_drsasa.hpp"

class DecoupledDSASA : public BaseDRSASA {
protected:
    int Imode;
    bool mtrx = true;
    
    void determine_Imode() {
        if (chain_sep.size() <= 1) {
            Imode = 2;  // Automatic molecular contact solver
        } else {
            Imode = 3;  // Automatic chain contact solver
        }
        
        // Check for protein chains
        set<string> proteinChains;
        if (chain_sep.size() <= 1) {
            for (const auto& atom: pdb) {
                if (atom.MOL_TYPE == "PROTEIN") {
                    proteinChains.insert(atom.CHAIN);
                } else {
                    proteinChains.clear();
                    break;
                }
            }
            
            if (proteinChains.size() == 2) {
                Imode = 3;  // Set to chain mode for two protein chains
            }
        }
    }

    void solve() override {
        determine_Imode();
        SolveInteractions(pdb, Imode);
        DecoupledSolver(pdb, rad.Points);
    }

    DRSASAResult process_results() override {
        DRSASAResult result;
        const size_t n_atoms = pdb.size();

        // Create basic arrays
        auto dsasa = py::array_t<double>(n_atoms);
        py::buffer_info dsasa_buf = dsasa.request();
        double* dsasa_ptr = static_cast<double*>(dsasa_buf.ptr);

        // Process contact data
        map<string, double> contact_by_type;
        map<string, double> contact_by_chain;
        map<pair<string, string>, double> chain_contacts;
        map<string, vector<size_t>> chain_atoms;
        
        double total_contact_area = 0.0;
        
        // First pass: collect basic data
        for(size_t i = 0; i < n_atoms; ++i) {
            const auto& atom = pdb[i];
            dsasa_ptr[i] = atom.dSASA;
            
            if(atom.dSASA > 0) {
                contact_by_type[atom.MOL_TYPE] += atom.dSASA;
                contact_by_chain[atom.CHAIN] += atom.dSASA;
                total_contact_area += atom.dSASA;
            }
            
            chain_atoms[atom.CHAIN].push_back(i);
        }

        // Process chain contacts
        if(Imode == 3) {  // Chain mode
            for(const auto& [chain1, atoms1] : chain_atoms) {
                for(const auto& [chain2, atoms2] : chain_atoms) {
                    if(chain1 < chain2) {  // Avoid double counting
                        double contact = 0.0;
                        for(size_t idx1 : atoms1) {
                            for(size_t idx2 : atoms2) {
                                for(const auto& interaction : pdb[idx1].interactions) {
                                    if(interaction.first == idx2) {
                                        contact += interaction.second;
                                    }
                                }
                            }
                        }
                        if(contact > 0) {
                            chain_contacts[{chain1, chain2}] = contact;
                        }
                    }
                }
            }
        }

        // Create interaction matrix if requested
        if(mtrx) {
            auto matrix = py::array_t<double>({n_atoms, n_atoms});
            py::buffer_info matrix_buf = matrix.request();
            double* matrix_ptr = static_cast<double*>(matrix_buf.ptr);

            // Fill interaction matrix
            for(size_t i = 0; i < n_atoms; ++i) {
                for(const auto& interaction : pdb[i].interactions) {
                    size_t j = interaction.first;
                    double value = interaction.second;
                    matrix_ptr[i * n_atoms + j] = value;
                    matrix_ptr[j * n_atoms + i] = value;  // Symmetric
                }
            }
            
            result.statistics["contact_matrix"] = matrix;
        }

        // Process molecular type contacts if in molecular mode
        map<pair<string, string>, double> mol_type_contacts;
        if(Imode == 2) {  // Molecular mode
            for(size_t i = 0; i < n_atoms; ++i) {
                const auto& atom1 = pdb[i];
                for(const auto& interaction : atom1.interactions) {
                    const auto& atom2 = pdb[interaction.first];
                    auto type_pair = make_pair(
                        min(atom1.MOL_TYPE, atom2.MOL_TYPE),
                        max(atom1.MOL_TYPE, atom2.MOL_TYPE)
                    );
                    mol_type_contacts[type_pair] += interaction.second;
                }
            }
        }

        // Pack results
        result.sasa = dsasa;
        result.metadata = py::dict(
            "contact_by_type"_a=contact_by_type,
            "contact_by_chain"_a=contact_by_chain,
            "chain_contacts"_a=chain_contacts,
            "molecular_contacts"_a=mol_type_contacts
        );
        
        result.statistics = py::dict(
            "total_contact_area"_a=total_contact_area,
            "mode"_a=Imode,
            "is_chain_mode"_a=(Imode == 3),
            "is_molecular_mode"_a=(Imode == 2)
        );

        return result;
    }

public:
    DecoupledDSASA() = default;
    
    void set_matrix_output(bool enable) { mtrx = enable; }
    
    // Utility methods
    bool is_chain_mode() const { return Imode == 3; }
    bool is_molecular_mode() const { return Imode == 2; }
    
    // Get specific contact areas
    py::dict get_chain_contacts() {
        if(!is_chain_mode()) {
            throw std::runtime_error("Chain contacts only available in chain mode");
        }
        return process_results().metadata["chain_contacts"].cast<py::dict>();
    }
    
    py::dict get_molecular_contacts() {
        if(!is_molecular_mode()) {
            throw std::runtime_error("Molecular contacts only available in molecular mode");
        }
        return process_results().metadata["molecular_contacts"].cast<py::dict>();
    }
};

PYBIND11_MODULE(drsasa, m) {
    m.doc() = "Python bindings for DR_SASA: Solvent Accessible Surface Area Calculator";
    
    // First bind the result structure
    py::class_<DRSASAResult>(m, "DRSASAResult")
        .def(py::init<>())
        .def_readonly("sasa", &DRSASAResult::sasa)
        .def_readonly("rel_sasa", &DRSASAResult::rel_sasa)
        .def_readonly("coordinates", &DRSASAResult::coordinates)
        .def_readonly("metadata", &DRSASAResult::metadata)
        .def_readonly("statistics", &DRSASAResult::statistics)
        .def("__repr__",
            [](const DRSASAResult& r) {
                return "DRSASAResult(total_atoms=" + 
                       std::to_string(r.sasa.size()) + ")";
            }
        );

    // Bind base class with common functionality
    py::class_<BaseDRSASA>(m, "BaseDRSASA")
        .def("set_probe_radius", &BaseDRSASA::set_probe_radius)
        .def("set_vdw_file", &BaseDRSASA::set_vdw_file)
        .def("set_chains", &BaseDRSASA::set_chains)
        .def("set_cl_mode", &BaseDRSASA::set_cl_mode)
        .def("enable_reorder", &BaseDRSASA::enable_reorder);

    // Mode 0: Simple SASA
    py::class_<SimpleSASA, BaseDRSASA>(m, "SimpleSASA")
        .def(py::init<>())
        .def("run", &SimpleSASA::run,
            py::arg("input_file"),
            "Calculate simple SASA for a structure")
        .def("__repr__",
            [](const SimpleSASA&) {
                return "SimpleSASA(mode=0)";
            }
        );

    // Mode 100: Relative SASA
    py::class_<RelativeSASA, SimpleSASA>(m, "RelativeSASA")
        .def(py::init<>())
        .def("run", &RelativeSASA::run,
            py::arg("input_file"),
            "Calculate relative SASA for a structure")
        .def("__repr__",
            [](const RelativeSASA&) {
                return "RelativeSASA(mode=100)";
            }
        );

    // Mode 1: Generic dSASA
    py::class_<GenericDSASA, BaseDRSASA>(m, "GenericDSASA")
        .def(py::init<>())
        .def("run", &GenericDSASA::run,
            py::arg("input_file"),
            "Calculate generic delta SASA between chains")
        .def("set_matrix_output", &GenericDSASA::set_matrix_output)
        .def("get_interaction_matrix", &GenericDSASA::get_interaction_matrix)
        .def("__repr__",
            [](const GenericDSASA&) {
                return "GenericDSASA(mode=1)";
            }
        );

    // Modes 2 & 3: Internal dSASA
    py::class_<InternalDSASA, BaseDRSASA>(m, "InternalDSASA")
        .def(py::init<bool>(), py::arg("residue_mode")=true)
        .def_static("create_residue_mode", &InternalDSASA::create_residue_mode)
        .def_static("create_atom_mode", &InternalDSASA::create_atom_mode)
        .def("run", &InternalDSASA::run,
            py::arg("input_file"),
            "Calculate internal delta SASA")
        .def("set_matrix_output", &InternalDSASA::set_matrix_output)
        .def("is_residue_mode", &InternalDSASA::is_residue_mode)
        .def("__repr__",
            [](const InternalDSASA& self) {
                return "InternalDSASA(mode=" + 
                       std::string(self.is_residue_mode() ? "2_residue" : "3_atom") + 
                       ")";
            }
        );

    // Mode 4: Decoupled dSASA
    py::class_<DecoupledDSASA, BaseDRSASA>(m, "DecoupledDSASA")
        .def(py::init<>())
        .def("run", &DecoupledDSASA::run,
            py::arg("input_file"),
            "Calculate decoupled delta SASA")
        .def("set_matrix_output", &DecoupledDSASA::set_matrix_output)
        .def("is_chain_mode", &DecoupledDSASA::is_chain_mode)
        .def("is_molecular_mode", &DecoupledDSASA::is_molecular_mode)
        .def("get_chain_contacts", &DecoupledDSASA::get_chain_contacts)
        .def("get_molecular_contacts", &DecoupledDSASA::get_molecular_contacts)
        .def("__repr__",
            [](const DecoupledDSASA&) {
                return "DecoupledDSASA(mode=4)";
            }
        );

    // Add convenience functions
    m.def("calculate_simple_sasa", [](const std::string& input_file) {
        SimpleSASA calculator;
        return calculator.run(input_file);
    }, "Quick calculation of simple SASA");

    m.def("calculate_delta_sasa", [](const std::string& input_file,
                                   const std::vector<std::vector<std::string>>& chains) {
        GenericDSASA calculator;
        calculator.set_chains(chains);
        return calculator.run(input_file);
    }, "Quick calculation of delta SASA between chains");

    // Version information
    m.attr("__version__") = "0.5.0";
    m.attr("__author__") = "Original: Ribeiro J., Ríos-Vera C., Melo F., Schüller A.";
}