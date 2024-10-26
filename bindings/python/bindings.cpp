// bindings/python/bindings.cpp
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <string>
#include <vector>
#include <map>
#include <stdexcept>
#include <cstdint>

#include "stdafx.h"
#include "atom_struct.h"
#include "PDBparser2.h"
#include "SetRadius.h"
#include "SurfaceSolverCL.h"

namespace py = pybind11;
using std::string;
using std::vector;
using std::map;
using std::runtime_error;

// External function declarations
extern void SimpleSolverCL(vector<atom_struct>& pdb, vector<float>& points, int cl_mode);
extern vector<atom_struct> PDBparser(const string& fname, const string& type_dict, bool keep_unknown);

class SimpleSASA {
public:
    SimpleSASA(float probe_radius = 1.4f) : probe_radius_(probe_radius) {
        vdw_radii_.GenPoints(0);  // Initialize with default points
    }

    py::dict calculate_sasa(const string& pdb_file) {
        try {
            // Parse PDB and calculate SASA
            auto atoms = PDBparser(pdb_file, "", true);
            if (atoms.empty()) {
                throw runtime_error("No atoms loaded from PDB file");
            }
            
            // Set radii and calculate SASA
            vdw_radii_.SetRadius(atoms, probe_radius_);
            SimpleSolverCL(atoms, vdw_radii_.Points, 0);  // Use CPU mode
            
            return create_results_dict(atoms);
        } catch (const std::exception& e) {
            throw runtime_error(string("SASA calculation failed: ") + e.what());
        }
    }

    float get_probe_radius() const { return probe_radius_; }
    void set_probe_radius(float radius) { probe_radius_ = radius; }

private:
    float probe_radius_;
    VDWcontainer vdw_radii_;

    py::dict create_results_dict(const vector<atom_struct>& atoms) const {
        py::dict results;
        
        // Prepare data arrays
        vector<double> atom_sasa;
        vector<string> atom_names;
        vector<string> residue_names;
        vector<int> residue_numbers;
        vector<string> chain_ids;
        vector<float> coordinates;
        
        double total_sasa = 0.0;
        
        // Populate arrays
        for (const auto& atom : atoms) {
            atom_sasa.push_back(atom.SASA);
            atom_names.push_back(atom.NAME);
            residue_names.push_back(atom.RESN);
            residue_numbers.push_back(atom.RESI);
            chain_ids.push_back(atom.CHAIN);
            
            coordinates.push_back(atom.COORDS[0]);
            coordinates.push_back(atom.COORDS[1]);
            coordinates.push_back(atom.COORDS[2]);
            
            total_sasa += atom.SASA;
        }
        
        // Create numpy arrays
        results["atom_sasa"] = py::array_t<double>(atom_sasa.size(), atom_sasa.data());
        results["atom_names"] = atom_names;
        results["residue_names"] = residue_names;
        results["residue_numbers"] = py::array_t<int>(residue_numbers.size(), residue_numbers.data());
        results["chain_ids"] = chain_ids;
        
        // Create coordinates as Nx3 array
        vector<ssize_t> shape = {static_cast<ssize_t>(atoms.size()), 3};
        vector<ssize_t> strides = {3 * sizeof(float), sizeof(float)};
        results["coordinates"] = py::array_t<float>(shape, strides, coordinates.data());
        
        results["total_sasa"] = total_sasa;
        
        return results;
    }
};

PYBIND11_MODULE(dr_sasa_py, m) {
    m.doc() = R"pbdoc(
        Python bindings for dr_sasa
        --------------------------
        A simple interface for calculating Solvent Accessible Surface Area (SASA)
    )pbdoc";

    py::class_<SimpleSASA>(m, "SimpleSASA")
        .def(py::init<float>(), py::arg("probe_radius") = 1.4f)
        .def("calculate_sasa", &SimpleSASA::calculate_sasa, py::arg("pdb_file"),
             "Calculate SASA for all atoms in the structure")
        .def_property("probe_radius", 
                     &SimpleSASA::get_probe_radius,
                     &SimpleSASA::set_probe_radius);

    m.attr("__version__") = "0.1.0";
}