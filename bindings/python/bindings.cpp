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
#include "SurfaceSolverOnTheFly.h"

namespace py = pybind11;
using namespace py::literals;  // Adds support for _a literal

using std::string;
using std::vector;
using std::map;
using std::set;
using std::runtime_error;

// External function declarations
extern void SimpleSolverCL(vector<atom_struct>& pdb, vector<float>& points, int cl_mode);
extern void Generic_Solver(vector<atom_struct>& pdb, vector<float>& points, vector<vector<string>> obj1, int mode, int cl_mode);
extern void DecoupledSolver(vector<atom_struct>& pdb, vector<float>& points);
extern void ChainSelector(vector<vector<string>>& selection, vector<atom_struct>& pdb);
extern void RelativeSASA(vector<atom_struct>& pdb);
extern void GeneratePairInteractionData(vector<atom_struct>& pdb);
extern void SolveInteractions(vector<atom_struct>& pdb, uint32 mode);
// Additional useful external functions we might need later
//extern void CalculateDNA_ProtInteractions(vector<atom_struct>& pdb, int mode);
//extern vector<atom_struct> ReorderPDB(const vector<atom_struct>& pdb);

static constexpr float DEFAULT_PROBE_RADIUS = 1.4f;  // Water probe in Angstroms

// Simple SASA (Mode 0)
class SimpleSASA {
private:
    float probe_radius_;
    int cl_mode_;
    VDWcontainer vdw_radii_;

public:
    SimpleSASA(float probe_radius = 1.4f, int compute_mode = 0) 
        : probe_radius_(probe_radius), cl_mode_(compute_mode) {
        vdw_radii_.GenPoints();
    }
    
    py::dict calculate(const string& pdb_file);
};

py::dict SimpleSASA::calculate(const string& pdb_file) {
    auto atoms = PDBparser(pdb_file, "", true);
    if (atoms.empty()) {
        throw std::runtime_error("No atoms loaded from PDB file");
    }

    vdw_radii_.SetRadius(atoms, probe_radius_);
    SolveInteractions(atoms, 0);  // Mode 0 specific
    SimpleSolverCL(atoms, vdw_radii_.Points, cl_mode_);
    
    return create_analysis_results(atoms, false);
}

// Generic dSASA (Mode 1)
class GenericSASA {
private:
    float probe_radius_;
    int cl_mode_;
    VDWcontainer vdw_radii_;

public:
    GenericSASA(float probe_radius = 1.4f, int compute_mode = 0) 
        : probe_radius_(probe_radius), cl_mode_(compute_mode) {
        vdw_radii_.GenPoints();
    }
    
    py::dict calculate(const string& pdb_file, vector<vector<string>>& chains);
};

py::dict GenericSASA::calculate(const string& pdb_file, vector<vector<string>>& chains) {  
    auto atoms = PDBparser(pdb_file, "", true);
    if (atoms.empty()) {
        throw std::runtime_error("No atoms loaded from PDB file");
    }

    vdw_radii_.SetRadius(atoms, probe_radius_);
    
    // Determine mode
    int Imode = (chains.size() <= 1) ? 4 : 1;  // 4=Auto, 1=Manual

    // Check for protein-only case
    if (chains.size() <= 1) {
        set<string> proteinChains;
        for (const auto& atom: atoms) {
            if (atom.MOL_TYPE == "PROTEIN") {
                proteinChains.insert(atom.CHAIN);
            } else {
                proteinChains.clear();
                break;
            }
        }
        if (proteinChains.size() == 2) {
            Imode = 5;  // Protein-protein mode
        }
    }

    ChainSelector(chains, atoms);
    Generic_Solver(atoms, vdw_radii_.Points, chains, Imode, cl_mode_);
    GeneratePairInteractionData(atoms);
    
    return create_analysis_results(atoms, true);
}

// Decoupled dSASA (Mode 4)
class DecoupledSASA {
private:
    float probe_radius_;
    int cl_mode_;
    VDWcontainer vdw_radii_;

public:
    DecoupledSASA(float probe_radius = 1.4f, int compute_mode = 0) 
        : probe_radius_(probe_radius), cl_mode_(compute_mode) {
        vdw_radii_.GenPoints();
    }
    
    py::dict calculate(const string& pdb_file, vector<vector<string>>& chains);
};

py::dict DecoupledSASA::calculate(const string& pdb_file, vector<vector<string>>& chains) {
    auto atoms = PDBparser(pdb_file, "", true);
    if (atoms.empty()) {
        throw std::runtime_error("No atoms loaded from PDB file");
    }

    vdw_radii_.SetRadius(atoms, probe_radius_);
    
    int Imode = (chains.size() <= 1) ? 2 : 3;  // 2=Molecular, 3=Chain

    if (!chains.empty()) {
        ChainSelector(chains, atoms);
    }

    SolveInteractions(atoms, Imode);
    DecoupledSolver(atoms, vdw_radii_.Points);
    
    return create_analysis_results(atoms, true);
}

// Relative SASA (Mode 100)
class RelSASA {
private:
    float probe_radius_;
    int cl_mode_;
    VDWcontainer vdw_radii_;

public:
    RelSASA(float probe_radius = 1.4f, int compute_mode = 0) 
        : probe_radius_(probe_radius), cl_mode_(compute_mode) {
        vdw_radii_.GenPoints();
    }
    
    py::dict calculate(const string& pdb_file);
};

py::dict RelSASA::calculate(const string& pdb_file) {
    auto atoms = PDBparser(pdb_file, "", true);
    if (atoms.empty()) {
        throw std::runtime_error("No atoms loaded from PDB file");
    }

    vdw_radii_.SetRadius(atoms, probe_radius_);
    SolveInteractions(atoms, 0);
    SimpleSolverCL(atoms, vdw_radii_.Points, cl_mode_);
    RelativeSASA(atoms);  // Additional step for relative SASA
    
    return create_analysis_results(atoms, false);
}

// Check SolverDataPorcessing for simplyfing this. 
py::dict create_analysis_results(const vector<atom_struct>& atoms, bool include_matrix = true) {
    py::dict results;
    const size_t n_atoms = atoms.size();

    // 1. Basic Arrays - Atom Information
    auto sasa = py::array_t<double>(n_atoms);
    auto dsasa = py::array_t<double>(n_atoms);
    auto rel_sasa = py::array_t<double>(n_atoms);
    vector<string> atom_names;
    vector<string> elements;
    vector<string> mol_types;
    vector<float> coordinates;
    vector<float> radii;
    vector<float> occupancies;
    vector<float> b_factors;
    vector<string> charges;
    vector<bool> is_hetatm;
    vector<int> atom_types;
    
    // Create buffers
    {
        py::buffer_info sasa_buf = sasa.request();
        py::buffer_info dsasa_buf = dsasa.request();
        py::buffer_info rel_buf = rel_sasa.request();
        
        double* sasa_ptr = static_cast<double*>(sasa_buf.ptr);
        double* dsasa_ptr = static_cast<double*>(dsasa_buf.ptr);
        double* rel_ptr = static_cast<double*>(rel_buf.ptr);
        
        for(size_t i = 0; i < n_atoms; ++i) {
            const auto& atom = atoms[i];
            
            // SASA information
            sasa_ptr[i] = atom.SASA;
            dsasa_ptr[i] = atom.EXT0;  // dSASA stored in EXT0
            rel_ptr[i] = atom.ACCESSIBLE_P;
            
            // Atom properties
            atom_names.push_back(atom.NAME);
            elements.push_back(atom.ELEMENT);
            mol_types.push_back(atom.MOL_TYPE);
            radii.push_back(atom.RADIUS);
            occupancies.push_back(atom.OCCUPANCY);
            b_factors.push_back(atom.TFACTOR);
            charges.push_back(atom.CHARGE);
            is_hetatm.push_back(atom.HETATM);
            atom_types.push_back(atom.ATOM_TYPE);
            
            // Coordinates
            coordinates.push_back(atom.COORDS[0]);
            coordinates.push_back(atom.COORDS[1]);
            coordinates.push_back(atom.COORDS[2]);
        }
    }

    // 2. Molecular Type Analysis
    map<string, double> sasa_by_mol_type;
    map<string, size_t> atoms_by_type;
    map<string, size_t> residues_by_type;
    map<pair<string, string>, double> mol_type_contacts;
    set<string> unique_mol_types;
    set<string> unique_residues;

    // 3. Interface Analysis
    map<string, vector<size_t>> interface_atoms;
    double total_interface_area = 0.0;
    map<string, double> interface_by_chain;
    double buried_surface_area = 0.0;

    // Process atoms
    for(size_t i = 0; i < n_atoms; ++i) {
        const auto& atom = atoms[i];
        string mol_type = atom.MOL_TYPE;
        string res_key = atom.CHAIN + "_" + atom.RESN + "_" + std::to_string(atom.RESI);
        
        // Molecular type stats
        sasa_by_mol_type[mol_type] += atom.SASA;
        atoms_by_type[mol_type]++;
        unique_mol_types.insert(mol_type);
        unique_residues.insert(res_key);
        residues_by_type[mol_type] = unique_residues.size();

        // Interface analysis
        if(atom.INT_NUM > 0) {  // If atom is at interface
            interface_atoms[atom.CHAIN].push_back(i);
            interface_by_chain[atom.CHAIN] += atom.EXT0;
            total_interface_area += atom.EXT0;
            buried_surface_area += (atom.AREA - atom.SASA);
        }

        // Process interactions and contact areas
        for(const auto& interaction : atom.INTERACTION_P) {
            const auto& partner = atoms[interaction];
            auto type_pair = make_pair(
                min(atom.MOL_TYPE, partner.MOL_TYPE),
                max(atom.MOL_TYPE, partner.MOL_TYPE)
            );
            mol_type_contacts[type_pair] += 1; //atom.CONTACT_AREA[interaction]; needs to be fixed
        }
    }

    // 4. Coordinates as Nx3 array
    vector<ssize_t> coord_shape = {static_cast<ssize_t>(n_atoms), 3};
    vector<ssize_t> coord_strides = {3 * sizeof(float), sizeof(float)};
    
    // Pack all results
    results["atom_info"] = py::dict(
        "names"_a=atom_names,
        "elements"_a=elements,
        "mol_types"_a=mol_types,
        "coordinates"_a=py::array_t<float>(coord_shape, coord_strides, coordinates.data()),
        "radii"_a=py::array_t<float>(radii.size(), radii.data()),
        "occupancies"_a=py::array_t<float>(occupancies.size(), occupancies.data()),
        "b_factors"_a=py::array_t<float>(b_factors.size(), b_factors.data()),
        "charges"_a=charges,
        "is_hetatm"_a=is_hetatm,
        "atom_types"_a=py::array_t<int>(atom_types.size(), atom_types.data())
    );

    results["sasa"] = py::dict(
        "values"_a=sasa,
        "delta"_a=dsasa,
        "relative"_a=rel_sasa,
        "by_mol_type"_a=sasa_by_mol_type
    );

    results["interface"] = py::dict(
        "total_area"_a=total_interface_area,
        "buried_surface_area"_a=buried_surface_area,
        "by_chain"_a=interface_by_chain,
        "atoms"_a=interface_atoms
    );

    results["molecular"] = py::dict(
        "atoms_by_type"_a=atoms_by_type,
        "residues_by_type"_a=residues_by_type,
        "type_contacts"_a=mol_type_contacts
    );

    // Include detailed interaction data if available
    if(!atoms.empty() && !atoms[0].INTERACTION_P.empty()) {
        map<pair<string, string>, vector<pair<size_t, size_t>>> interactions_by_type;
        vector<double> interaction_energies;
        
        for(size_t i = 0; i < n_atoms; ++i) {
            const auto& atom = atoms[i];
            for(const auto& j : atom.INTERACTION_P) {
                const auto& partner = atoms[j];
                auto type_pair = make_pair(atom.MOL_TYPE, partner.MOL_TYPE);
                interactions_by_type[type_pair].push_back({i, j});
                interaction_energies.push_back(atom.ENERGY);
            }
        }
        
        results["interactions"] = py::dict(
            "by_type"_a=interactions_by_type,
            "energies"_a=py::array_t<double>(interaction_energies.size(), interaction_energies.data())
        );
    }

    return results;
}

// Module definition
PYBIND11_MODULE(dr_sasa_py, m) {
    m.doc() = "Python bindings for DR_SASA: Solvent Accessible Surface Area Calculator";
    
    // Simple SASA (Mode 0)
    py::class_<SimpleSASA>(m, "SimpleSASA")
        .def(py::init<float, int>(),
             py::arg("probe_radius") = 1.4f,
             py::arg("compute_mode") = 0)
        .def("calculate", &SimpleSASA::calculate,
             py::arg("pdb_file"),
             "Calculate simple SASA using mode 0.\n"
             "Returns dict with SASA values and basic statistics.");

    // Generic dSASA (Mode 1)
    py::class_<GenericSASA>(m, "GenericSASA")
        .def(py::init<float, int>(),
             py::arg("probe_radius") = 1.4f,
             py::arg("compute_mode") = 0)
        .def("calculate", &GenericSASA::calculate,
             py::arg("pdb_file"),
             py::arg("chains"),
             "Calculate delta SASA between chain groups using mode 1.\n"
             "Automatically determines interaction mode (1=Manual, 4=Auto, 5=Protein-protein).\n"
             "Returns dict with SASA values, interactions, and matrices.");

    // Decoupled dSASA (Mode 4)
    py::class_<DecoupledSASA>(m, "DecoupledSASA")
        .def(py::init<float, int>(),
             py::arg("probe_radius") = 1.4f,
             py::arg("compute_mode") = 0)
        .def("calculate", &DecoupledSASA::calculate,
             py::arg("pdb_file"),
             py::arg("chains"),
             "Calculate decoupled delta SASA using mode 4.\n"
             "Automatically determines mode (2=Molecular, 3=Chain).\n"
             "Returns dict with SASA values and contact matrices.");

    // Relative SASA (Mode 100)
    py::class_<RelSASA>(m, "RelativeSASA")
        .def(py::init<float, int>(),
             py::arg("probe_radius") = 1.4f,
             py::arg("compute_mode") = 0)
        .def("calculate", &RelSASA::calculate,
             py::arg("pdb_file"),
             "Calculate relative SASA using mode 100.\n"
             "Returns dict with SASA values and relative accessibility.");

    // Convenience functions
    m.def("calculate_simple_sasa", [](const string& pdb_file, float probe_radius = 1.4f) {
        SimpleSASA calculator(probe_radius);
        return calculator.calculate(pdb_file);
    }, py::arg("pdb_file"), py::arg("probe_radius") = 1.4f,
    "Quick calculation of simple SASA (Mode 0)");

    m.def("calculate_delta_sasa", [](const string& pdb_file, 
                                    vector<vector<string>>& chains,  // Non-const chains!
                                    float probe_radius = 1.4f) {
        GenericSASA calculator(probe_radius);
        return calculator.calculate(pdb_file, chains);
    }, py::arg("pdb_file"), py::arg("chains"), py::arg("probe_radius") = 1.4f,
    "Quick calculation of delta SASA between chains (Mode 1)");

    // Module info
    m.attr("__version__") = "0.5.0";
    m.attr("__author__") = "Original: Ribeiro J., Ríos-Vera C., Melo F., Schüller A.";
    // Constants
    m.attr("DEFAULT_PROBE_RADIUS") = 1.4f;
}