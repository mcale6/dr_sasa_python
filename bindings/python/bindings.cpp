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

#include "stdafx.h"
#include "atom_struct.h"
#include "PDBparser2.h"
#include "SetRadius.h"
#include "SurfaceSolverCL.h"

namespace py = pybind11;

using std::string;
using std::vector;
using std::map;
using std::set;
using std::runtime_error;

// External function declarations
extern void SimpleSolverCL(vector<atom_struct>& pdb, vector<float>& points, int cl_mode);
extern void SimpleSolverCL(vector<atom_struct>& pdb, vector<float>& points, int cl_mode);
extern void Generic_Solver(vector<atom_struct>& pdb, vector<float>& points, vector<vector<string>> obj1, int mode, int cl_mode);
extern void DecoupledSolver(vector<atom_struct>& pdb, vector<float>& points);
extern void ChainSelector(vector<vector<string>>& chain_sep, vector<atom_struct>& pdb);
extern void RelativeSASA(vector<atom_struct>& pdb);
extern void GeneratePairInteractionData(vector<atom_struct>& pdb);

// Additional useful external functions we might need
extern void SolveInteractions(vector<atom_struct>& pdb, int mode);
extern void CalculateDNA_ProtInteractions(vector<atom_struct>& pdb, int mode);
extern vector<atom_struct> ReorderPDB(const vector<atom_struct>& pdb);

static constexpr float DEFAULT_PROBE_RADIUS = 1.4f;  // Water probe in Angstroms

class BaseSASA {
protected:
    float probe_radius_;
    int cl_mode_;
    VDWcontainer vdw_radii_;

    // Constructor in base class
    BaseSASA(float probe_radius, int compute_mode) 
        : probe_radius_(probe_radius), 
          cl_mode_(compute_mode) {
        vdw_radii_.GenPoints(0);
    }

    // Common results creation
    py::dict create_unified_results(const vector<atom_struct>& atoms) const;
};

// Simple SASA calculator
class SimpleSASA : public BaseSASA {
public:
    // Properly call base class constructor
    SimpleSASA(float probe_radius = DEFAULT_PROBE_RADIUS, int compute_mode = 0) 
        : BaseSASA(probe_radius, compute_mode) {}

    py::dict calculate(const string& pdb_file);
};

class GenericSASA : public BaseSASA {
public:
    GenericSASA(float probe_radius = DEFAULT_PROBE_RADIUS, int compute_mode = 0) 
        : BaseSASA(probe_radius, compute_mode) {}

    py::dict calculate(const string& pdb_file, 
                      const vector<vector<string>>& chains,
                      int mode);
};

class DecoupledSASA : public BaseSASA {
public:
    DecoupledSASA(float probe_radius = DEFAULT_PROBE_RADIUS, int compute_mode = 0) 
        : BaseSASA(probe_radius, compute_mode) {}

    py::dict calculate(const string& pdb_file);
};

// Method implementations
py::dict SimpleSASA::calculate(const string& pdb_file) {
    try {
        auto atoms = PDBparser(pdb_file, "", true);
        if (atoms.empty()) {
            throw runtime_error("No atoms loaded from PDB file");
        }
        
        vdw_radii_.SetRadius(atoms, probe_radius_);
        SimpleSolverCL(atoms, vdw_radii_.Points, cl_mode_);
        
        return create_unified_results(atoms);
    } catch (const std::exception& e) {
        throw runtime_error(string("SASA calculation failed: ") + e.what());
    }
}

py::dict GenericSASA::calculate(const string& pdb_file, 
                              const vector<vector<string>>& chains,
                              int mode) {
    try {
        auto atoms = PDBparser(pdb_file, "", true);
        if (atoms.empty()) {
            throw runtime_error("No atoms loaded from PDB file");
        }
        
        vdw_radii_.SetRadius(atoms, probe_radius_);
        Generic_Solver(atoms, vdw_radii_.Points, chains, mode, cl_mode_);
        
        return create_unified_results(atoms);
    } catch (const std::exception& e) {
        throw runtime_error(string("Generic SASA calculation failed: ") + e.what());
    }
}

py::dict DecoupledSASA::calculate(const string& pdb_file) {
    try {
        auto atoms = PDBparser(pdb_file, "", true);
        if (atoms.empty()) {
            throw runtime_error("No atoms loaded from PDB file");
        }
        
        vdw_radii_.SetRadius(atoms, probe_radius_);
        DecoupledSolver(atoms, vdw_radii_.Points);
        
        return create_unified_results(atoms);
    } catch (const std::exception& e) {
        throw runtime_error(string("Decoupled SASA calculation failed: ") + e.what());
    }
}

py::dict BaseSASA::create_unified_results(const vector<atom_struct>& atoms) const {
    py::dict results;
    
    // Basic properties
    vector<double> atom_sasa;
    vector<double> relative_sasa;
    vector<string> atom_names;
    vector<string> residue_names;
    vector<int> residue_numbers;
    vector<string> chain_ids;
    vector<string> elements;
    vector<string> mol_types;
    vector<float> coordinates;
    
    // Accumulate results
    map<string, double> sasa_by_type;
    double total_sasa = 0.0;
    
    for (const auto& atom : atoms) {
        // Basic properties
        atom_sasa.push_back(atom.SASA);
        atom_names.push_back(atom.NAME);
        residue_names.push_back(atom.RESN);
        residue_numbers.push_back(atom.RESI);
        chain_ids.push_back(atom.CHAIN);
        elements.push_back(atom.ELEMENT);
        mol_types.push_back(atom.MOL_TYPE);
        
        coordinates.push_back(atom.COORDS[0]);
        coordinates.push_back(atom.COORDS[1]);
        coordinates.push_back(atom.COORDS[2]);
        
        // Accumulate totals
        total_sasa += atom.SASA;
        sasa_by_type[atom.MOL_TYPE] += atom.SASA;
    }
    
    // Pack results
    results["atom_sasa"] = py::array_t<double>(atom_sasa.size(), atom_sasa.data());
    results["atom_names"] = atom_names;
    results["residue_names"] = residue_names;
    results["residue_numbers"] = py::array_t<int>(residue_numbers.size(), residue_numbers.data());
    results["chain_ids"] = chain_ids;
    results["elements"] = elements;
    results["mol_types"] = mol_types;
    
    // Coordinates as Nx3 array
    vector<ssize_t> shape = {static_cast<ssize_t>(atoms.size()), 3};
    vector<ssize_t> strides = {3 * sizeof(float), sizeof(float)};
    results["coordinates"] = py::array_t<float>(shape, strides, coordinates.data());
    
    // Summary statistics
    results["total_sasa"] = total_sasa;
    results["sasa_by_type"] = sasa_by_type;
    
    return results;
};

// Module definition
PYBIND11_MODULE(dr_sasa_py, m) {
    m.doc() = "Python bindings for dr_sasa - SASA calculation library";

    py::class_<SimpleSASA>(m, "SimpleSASA")
        .def(py::init<float, int>(),
             py::arg("probe_radius") = 1.4f,
             py::arg("compute_mode") = 0)
        .def("calculate", &SimpleSASA::calculate,
             py::arg("pdb_file"),
             "Calculate SASA using simple solver");

    py::class_<GenericSASA>(m, "GenericSASA")
        .def(py::init<float, int>(),
             py::arg("probe_radius") = 1.4f,
             py::arg("compute_mode") = 0)
        .def("calculate", &GenericSASA::calculate,
             py::arg("pdb_file"),
             py::arg("chains"),
             py::arg("mode") = 1,
             "Calculate SASA with chain-specific analysis");

    py::class_<DecoupledSASA>(m, "DecoupledSASA")
        .def(py::init<float, int>(),
             py::arg("probe_radius") = 1.4f,
             py::arg("compute_mode") = 0)
        .def("calculate", &DecoupledSASA::calculate,
             py::arg("pdb_file"),
             "Calculate SASA using decoupled solver");

    m.attr("__version__") = "0.1.0";
}