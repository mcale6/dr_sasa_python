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
#include "histogram.h"
#include "SolverDataProcessing.h"

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
// Matrix generation external functions
extern void GenerateInterBSAMatrix(vector<atom_struct>& pdb,
                                 map<vector<string>, vector<float>>& matrixIJatom,
                                 map<vector<string>, vector<float>>& matrixIJres,
                                 map<vector<string>, vector<string>>& COLatom,
                                 map<vector<string>, vector<string>>& COLres, 
                                 map<vector<string>, vector<string>>& ROWatom,
                                 map<vector<string>, vector<string>>& ROWres,
                                 map<vector<string>, vector<uint32>>& COLatomtype,
                                 map<vector<string>, vector<uint32>>& ROWatomtype);

extern void GenerateIntraBSAMatrix(vector<atom_struct>& pdb,
                                 vector<float>& matrixIJatom,
                                 vector<float>& matrixIJres,
                                 vector<string>& COLatom,
                                 vector<string>& COLres,
                                 vector<string>& ROWatom,
                                 vector<string>& ROWres,
                                 vector<uint32>& COLatomtype,
                                 vector<uint32>& ROWatomtype);
// Additional useful external functions we might need later
// GenerateInterBSAMatrix 
// GenerateIntraBSAMatrix
// CalculateDNA_ProtInteractions DNA /Prtoein?
// GeneratePDistribution
// AddRawAtomxNucData

// - SimpleSASA      (Mode 0)  # Basic SASA calculation
// - GenericSASA     (Mode 1)  # Chain/molecular type interactions
// - DecoupledSASA   (Mode 4)  # Decoupled interface analysis
// - RelativeSASA    (Mode 100) # Relative accessibility

// "Missing"
// - Mode 2 (Residue dSASA) # returned in the dict handled by GenericSASA/DecoupledSASA
// - Mode 3 (Atom dSASA)    # returned in the dicth andled by GenericSASA/DecoupledSASA

static constexpr float DEFAULT_PROBE_RADIUS = 1.4f;  // Water probe in Angstroms

// Check SolverDataPorcessing and atom.struct for more information.  
py::dict create_analysis_results(vector<atom_struct>& atoms, bool include_matrix = true) {
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
    vector<ssize_t> coord_shape = {static_cast<ssize_t>(n_atoms), 3};
    vector<ssize_t> coord_strides = {3 * sizeof(float), sizeof(float)};
    // Pack atom information
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

    // Add SASA results
    map<string, double> sasa_by_mol_type;
    for(const auto& atom : atoms) {
        sasa_by_mol_type[atom.MOL_TYPE] += atom.SASA;
    }

    results["sasa"] = py::dict(
        "values"_a=sasa,
        "delta"_a=dsasa,
        "relative"_a=rel_sasa,
        "by_mol_type"_a=sasa_by_mol_type
    );

    // Interface Analysis
    map<string, vector<size_t>> interface_atoms;
    double total_interface_area = 0.0;
    map<string, double> interface_by_chain;
    double buried_surface_area = 0.0;

    for(size_t i = 0; i < n_atoms; ++i) {
        const auto& atom = atoms[i];
        if(atom.INT_NUM > 0) {  // If atom is at interface
            interface_atoms[atom.CHAIN].push_back(i);
            interface_by_chain[atom.CHAIN] += atom.EXT0;
            total_interface_area += atom.EXT0;
            buried_surface_area += (atom.AREA - atom.SASA);
        }
    }

    results["interface"] = py::dict(
        "total_area"_a=total_interface_area,
        "buried_surface_area"_a=buried_surface_area,
        "by_chain"_a=interface_by_chain,
        "atoms"_a=interface_atoms
    );

    // Add matrix results if requested
    if(include_matrix && !atoms.empty()) {
        // Matrix containers for inter-molecular results
        map<vector<string>, vector<float>> matrixIJatom;
        map<vector<string>, vector<float>> matrixIJres;
        map<vector<string>, vector<string>> COLatom;
        map<vector<string>, vector<string>> COLres;
        map<vector<string>, vector<string>> ROWatom;
        map<vector<string>, vector<string>> ROWres;
        map<vector<string>, vector<uint32>> COLatomtype;
        map<vector<string>, vector<uint32>> ROWatomtype;

        // Generate inter-molecular matrices
        GenerateInterBSAMatrix(atoms, matrixIJatom, matrixIJres, 
                             COLatom, COLres, ROWatom, ROWres,
                             COLatomtype, ROWatomtype);

        // Containers for intra-molecular results
        vector<float> intraMatrixAtom;
        vector<float> intraMatrixRes;
        vector<string> intraColAtom;
        vector<string> intraColRes;
        vector<string> intraRowAtom;
        vector<string> intraRowRes;
        vector<uint32> intraColAtomType;
        vector<uint32> intraRowAtomType;

        // Generate intra-molecular matrices
        GenerateIntraBSAMatrix(atoms, intraMatrixAtom, intraMatrixRes,
                             intraColAtom, intraColRes, intraRowAtom, intraRowRes,
                             intraColAtomType, intraRowAtomType);

        // Add matrices to results
        results["matrices"] = py::dict(
            "inter_molecular"_a=py::dict(
                "atomic"_a=matrixIJatom,
                "residue"_a=matrixIJres,
                "col_atoms"_a=COLatom,
                "row_atoms"_a=ROWatom,
                "col_types"_a=COLatomtype,
                "row_types"_a=ROWatomtype
            ),
            "intra_molecular"_a=py::dict(
                "atomic"_a=intraMatrixAtom,
                "residue"_a=intraMatrixRes,
                "col_atoms"_a=intraColAtom,
                "row_atoms"_a=intraRowAtom,
                "col_types"_a=intraColAtomType,
                "row_types"_a=intraRowAtomType
            )
        );
    }
    return results;
}

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
    
    py::dict calculate(const string& pdb_file, bool include_matrix = true);
};

py::dict SimpleSASA::calculate(const string& pdb_file, bool include_matrix) {
    auto atoms = PDBparser(pdb_file, "", true);
    if (atoms.empty()) {
        throw std::runtime_error("No atoms loaded from PDB file");
    }

    vdw_radii_.SetRadius(atoms, probe_radius_);
    SolveInteractions(atoms, 0);  // Mode 0 specific
    SimpleSolverCL(atoms, vdw_radii_.Points, cl_mode_);
    
    return create_analysis_results(atoms, include_matrix);
}

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
    
    py::dict calculate(const string& pdb_file, vector<vector<string>>& chains, bool include_matrix = true);
};

py::dict GenericSASA::calculate(const string& pdb_file, vector<vector<string>>& chains, bool include_matrix) {
    auto atoms = PDBparser(pdb_file, "", true);
    if (atoms.empty()) {
        throw std::runtime_error("No atoms loaded from PDB file");
    }

    vdw_radii_.SetRadius(atoms, probe_radius_);
    
    // Convert Python tuple/list into proper C++ vector format
    vector<vector<string>> cpp_chains;
    for (const auto& chain_group : chains) {
        vector<string> group;
        for (const auto& chain : chain_group) {
            group.push_back(string(chain));
        }
        cpp_chains.push_back(group);
    }
    
    int Imode = (cpp_chains.size() <= 1) ? 4 : 1;  // 4=Auto, 1=Manual

    // Check for protein-only case
    if (cpp_chains.size() <= 1) {
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

    ChainSelector(cpp_chains, atoms);
    Generic_Solver(atoms, vdw_radii_.Points, cpp_chains, Imode, cl_mode_);
    GeneratePairInteractionData(atoms);
    
    return create_analysis_results(atoms, include_matrix);
}

// DecoupledSASA class
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
    
    py::dict calculate(const string& pdb_file, vector<vector<string>>& chains, bool include_matrix = true);
};

py::dict DecoupledSASA::calculate(const string& pdb_file, vector<vector<string>>& chains, bool include_matrix) {
    auto atoms = PDBparser(pdb_file, "", true);
    if (atoms.empty()) {
        throw std::runtime_error("No atoms loaded from PDB file");
    }

    vdw_radii_.SetRadius(atoms, probe_radius_);
    
    // Convert Python tuple/list into proper C++ vector format
    vector<vector<string>> cpp_chains;
    for (const auto& chain_group : chains) {
        vector<string> group;
        for (const auto& chain : chain_group) {
            group.push_back(string(chain));
        }
        cpp_chains.push_back(group);
    }
    
    int Imode = (cpp_chains.size() <= 1) ? 2 : 3;  // 2=Molecular, 3=Chain

    if (!cpp_chains.empty()) {
        ChainSelector(cpp_chains, atoms);
    }

    SolveInteractions(atoms, Imode);
    DecoupledSolver(atoms, vdw_radii_.Points);
    
    return create_analysis_results(atoms, include_matrix);
}


// RelSASA class
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
    
    py::dict calculate(const string& pdb_file, bool include_matrix = false);
};

py::dict RelSASA::calculate(const string& pdb_file, bool include_matrix) {
    auto atoms = PDBparser(pdb_file, "", true);
    if (atoms.empty()) {
        throw std::runtime_error("No atoms loaded from PDB file");
    }

    vdw_radii_.SetRadius(atoms, probe_radius_);
    SolveInteractions(atoms, 0);
    SimpleSolverCL(atoms, vdw_radii_.Points, cl_mode_);
    RelativeSASA(atoms);  // Additional step for relative SASA
    
    return create_analysis_results(atoms, include_matrix);
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
             py::arg("include_matrix") = true,
             "Calculate simple SASA using mode 0.\n"
             "Args:\n"
             "    pdb_file: Path to PDB file\n"
             "    include_matrix: Whether to include interaction matrices (default: True)\n"
             "Returns dict with SASA values, basic statistics, and optional matrices.");

    // Generic dSASA (Mode 1)
    py::class_<GenericSASA>(m, "GenericSASA")
        .def(py::init<float, int>(),
             py::arg("probe_radius") = 1.4f,
             py::arg("compute_mode") = 0)
        .def("calculate", &GenericSASA::calculate,
             py::arg("pdb_file"),
             py::arg("chains") = vector<vector<string>>(),  // Default empty vector
             py::arg("include_matrix") = true,
             "Calculate delta SASA between chain groups using mode 1.\n"
             "Args:\n"
             "    pdb_file: Path to PDB file\n"
             "    chains: List of chain groups\n"
             "    include_matrix: Whether to include interaction matrices (default: True)\n"
             "Automatically determines interaction mode (1=Manual, 4=Auto, 5=Protein-protein).\n"
             "Returns dict with SASA values, interactions, and optional matrices.");

    // Decoupled dSASA (Mode 4)
    py::class_<DecoupledSASA>(m, "DecoupledSASA")
        .def(py::init<float, int>(),
             py::arg("probe_radius") = 1.4f,
             py::arg("compute_mode") = 0)
        .def("calculate", &DecoupledSASA::calculate,
             py::arg("pdb_file"),
             py::arg("chains") = vector<vector<string>>(),  // Default empty vector
             py::arg("include_matrix") = true,
             "Calculate decoupled delta SASA using mode 4.\n"
             "Args:\n"
             "    pdb_file: Path to PDB file\n"
             "    chains: List of chain groups\n"
             "    include_matrix: Whether to include interaction matrices (default: True)\n"
             "Automatically determines mode (2=Molecular, 3=Chain).\n"
             "Returns dict with SASA values, contact matrices, and optional matrices.");

    // Relative SASA (Mode 100)
    py::class_<RelSASA>(m, "RelativeSASA")
        .def(py::init<float, int>(),
             py::arg("probe_radius") = 1.4f,
             py::arg("compute_mode") = 0)
        .def("calculate", &RelSASA::calculate,
             py::arg("pdb_file"),
             py::arg("include_matrix") = true,
             "Calculate relative SASA using mode 100.\n"
             "Args:\n"
             "    pdb_file: Path to PDB file\n"
             "    include_matrix: Whether to include interaction matrices (default: True)\n"
             "Returns dict with SASA values, relative accessibility, and optional matrices.");

    // Convenience functions
    m.def("calculate_simple_sasa", 
        [](const string& pdb_file, float probe_radius = 1.4f, bool include_matrix = true) {
            SimpleSASA calculator(probe_radius);
            return calculator.calculate(pdb_file, include_matrix);
        }, 
        py::arg("pdb_file"), 
        py::arg("probe_radius") = 1.4f,
        py::arg("include_matrix") = true,
        "Quick calculation of simple SASA (Mode 0)");

    m.def("calculate_delta_sasa", 
        [](const string& pdb_file, 
           vector<vector<string>>& chains,  
           float probe_radius = 1.4f,
           bool include_matrix = true) {
            GenericSASA calculator(probe_radius);
            return calculator.calculate(pdb_file, chains, include_matrix);
        }, 
        py::arg("pdb_file"), 
        py::arg("chains"), 
        py::arg("probe_radius") = 1.4f,
        py::arg("include_matrix") = true,
        "Quick calculation of delta SASA between chains (Mode 1)");

    // Module info
    m.attr("__version__") = "0.5.0";
    m.attr("__author__") = "Original: Ribeiro J., Ríos-Vera C., Melo F., Schüller A.";
    // Constants
    m.attr("DEFAULT_PROBE_RADIUS") = 1.4f;
}