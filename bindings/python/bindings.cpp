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
// GenerateInterBSAMatrix  for overlaps with CalculateDNA_ProtInteractions DNA /Prtoein?
// GenerateIntraBSAMatrix for overlaps with CalculateDNA_ProtInteractions DNA /Prtoein?
// PrintDNA_ProtResults for overlaps with CalculateDNA_ProtInteractions DNA /Prtoein?
// GeneratePDistribution 

// - SimpleSASA      (Mode 0)  # Basic SASA calculation
// - GenericSASA     (Mode 1)  # Chain/molecular type interactions
// - DecoupledSASA   (Mode 4)  # Decoupled interface analysis
// - RelativeSASA    (Mode 100) # Relative accessibility this is wrong in the SolverDataProcessing per atom not per residue!
 
// "Missing"
// - Mode 2 (Residue dSASA) # returned in the dict handled by GenericSASA/DecoupledSASA
// - Mode 3 (Atom dSASA)    # returned in the dict andled by GenericSASA/DecoupledSASA

static constexpr float DEFAULT_PROBE_RADIUS = 1.4f;  // Water probe in Angstroms

// Check SolverDataPorcessing and atom.struct for more information.  
py::dict create_analysis_results(vector<atom_struct>& atoms, bool include_matrix = true) {
    py::dict results;
    const size_t n_atoms = atoms.size();

    // Basic atom info arrays
    vector<string> structures, names, residues, chains, elements;
    vector<string> struct_types, mol_types, altlocs, icodes, dtypes;
    vector<int> residue_nums, atom_ids;
    vector<bool> is_hetatm, is_active, is_terminal;
    vector<int> atom_types, atom_types_40, int_nums, ef_int_nums;
    vector<int> polarity;

    // Geometric properties
    vector<float> coordinates, radii, radii2, vdw;
    vector<float> occupancies, b_factors;
    vector<string> charges;
    
    // Surface areas and volumes
    auto sasa = py::array_t<double>(n_atoms);
    auto area = py::array_t<double>(n_atoms);
    auto ext0 = py::array_t<double>(n_atoms);
    auto ext1 = py::array_t<double>(n_atoms);
    auto accessible = py::array_t<double>(n_atoms);
    
    // Energy terms
    vector<float> energies, dd_energies, as_energies;
    vector<float> cs_energies, sdd_energies, bsa_energies;

    // Interaction data 
    vector<vector<uint32_t>> interaction_atoms;
    vector<vector<uint32_t>> interaction_sasa_atoms;
    vector<map<uint32_t, float>> contact_areas;
    vector<vector<float>> distances;
    vector<vector<float>> dd_distances;
    vector<bool> dd_interacting;
        
    // Buried area data
    vector<vector<vector<uint32_t>>> buried_by_atoms;
    vector<vector<uint32_t>> buried_by_atoms_valid;
    vector<vector<float>> buried_areas;
    
    // Bonded atoms
    vector<vector<uint32_t>> bonded_atoms;
    
    // Overlap data
    vector<vector<vector<uint32_t>>> overlap_tables;
    vector<vector<float>> overlap_areas;
    vector<vector<float>> normalized_overlap_areas;

    {
        // Setup buffer pointers
        py::buffer_info sasa_buf = sasa.request();
        py::buffer_info area_buf = area.request();
        py::buffer_info ext0_buf = ext0.request();
        py::buffer_info ext1_buf = ext1.request();
        py::buffer_info acc_buf = accessible.request();
        
        double* sasa_ptr = static_cast<double*>(sasa_buf.ptr);
        double* area_ptr = static_cast<double*>(area_buf.ptr);
        double* ext0_ptr = static_cast<double*>(ext0_buf.ptr);
        double* ext1_ptr = static_cast<double*>(ext1_buf.ptr);
        double* acc_ptr = static_cast<double*>(acc_buf.ptr);
        
        for(size_t i = 0; i < n_atoms; ++i) {
            const auto& atom = atoms[i];
            
            // Basic properties
            structures.push_back(atom.STRUCTURE);
            atom_ids.push_back(i);
            names.push_back(atom.NAME);
            residues.push_back(atom.RESN);
            chains.push_back(atom.CHAIN);
            elements.push_back(atom.ELEMENT);
            struct_types.push_back(atom.STRUCT_TYPE);
            mol_types.push_back(atom.MOL_TYPE);
            altlocs.push_back(atom.ALTLOC);
            icodes.push_back(atom.iCODE);
            residue_nums.push_back(atom.RESI);
            
            // Flags
            is_hetatm.push_back(atom.HETATM);
            is_active.push_back(atom.ACTIVE);
            is_terminal.push_back(atom.TERMINAL);
            
            // Types
            atom_types.push_back(atom.ATOM_TYPE);
            atom_types_40.push_back(atom.ATOM_TYPE_40);
            int_nums.push_back(atom.INT_NUM);
            ef_int_nums.push_back(atom.EF_INT_NUM);
            dtypes.push_back(atom.DTYPE);
            polarity.push_back(atom.POLAR);
            
            // Geometric properties
            coordinates.insert(coordinates.end(), atom.COORDS.begin(), atom.COORDS.end());
            radii.push_back(atom.RADIUS);
            radii2.push_back(atom.RADIUS2);
            vdw.push_back(atom.VDW);
            occupancies.push_back(atom.OCCUPANCY);
            b_factors.push_back(atom.TFACTOR);
            charges.push_back(atom.CHARGE);
            
            // Surface areas
            sasa_ptr[i] = atom.SASA;
            area_ptr[i] = atom.AREA;
            ext0_ptr[i] = atom.EXT0;
            ext1_ptr[i] = atom.EXT1;
            acc_ptr[i] = atom.ACCESSIBLE_P;
            
            // Energy terms
            energies.push_back(atom.ENERGY);
            dd_energies.push_back(atom.DD_ENERGY);
            as_energies.push_back(atom.AS_ENERGY);
            cs_energies.push_back(atom.CS_ENERGY);
            sdd_energies.push_back(atom.sDD_ENERGY);
            bsa_energies.push_back(atom.BSA_ENERGY);
            
            // Interaction data
            interaction_atoms.push_back(atom.INTERACTION_P);
            interaction_sasa_atoms.push_back(atom.INTERACTION_SASA_P);
            contact_areas.push_back(atom.CONTACT_AREA);
            //distances.push_back(atom.DISTANCES);
            //dd_distances.push_back(atom.DD_DISTANCES);
            //dd_interacting.push_back(atom.DD_INTERACTING);
            
              
            // Buried area data
            buried_by_atoms.push_back(atom.AREA_BURIED_BY_ATOM_vector);
            buried_by_atoms_valid.push_back(atom.AREA_BURIED_BY_ATOM_vector_valid);
            buried_areas.push_back(atom.AREA_BURIED_BY_ATOM_area);
            
            // Bonded atoms
            bonded_atoms.push_back(atom.BONDED);
            
            // Overlap data
            overlap_tables.push_back(atom.ov_table);
            overlap_areas.push_back(atom.ov_table_area);
            normalized_overlap_areas.push_back(atom.ov_norm_area);
        }
    }

    // Pack coordinates with proper shape
    vector<ssize_t> coord_shape = {static_cast<ssize_t>(n_atoms), 3};
    vector<ssize_t> coord_strides = {3 * sizeof(float), sizeof(float)};

    // Construct comprehensive results dictionary
    results["atoms"] = py::dict(
        "structure"_a=structures,
        "ids"_a=atom_ids,
        "names"_a=names,
        "residues"_a=residues,
        "chains"_a=chains,
        "elements"_a=elements,
        "struct_types"_a=struct_types,
        "mol_types"_a=mol_types,
        "altlocs"_a=altlocs,
        "icodes"_a=icodes,
        "residue_numbers"_a=residue_nums,
        "is_hetatm"_a=is_hetatm,
        "is_active"_a=is_active,
        "is_terminal"_a=is_terminal,
        "atom_types"_a=atom_types,
        "atom_types_40"_a=atom_types_40,
        "int_nums"_a=int_nums,
        "ef_int_nums"_a=ef_int_nums,
        "dtypes"_a=dtypes,
        "polarity"_a=polarity,
        "coordinates"_a=py::array_t<float>(coord_shape, coord_strides, coordinates.data()),
        "radii"_a=radii,
        "radii2"_a=radii2,
        "vdw"_a=vdw,
        "occupancies"_a=occupancies,
        "b_factors"_a=b_factors,
        "charges"_a=charges,
        "bonded"_a=bonded_atoms
    );

    results["surface"] = py::dict(
        "area"_a=area,
        "sasa"_a=sasa,
        "buried_area"_a=ext0,
        "fast_dsasa"_a=ext1,
        "accessible"_a=accessible
    );

    results["energy"] = py::dict(
        "total"_a=energies,
        "dd"_a=dd_energies,
        "as"_a=as_energies,
        "cs"_a=cs_energies,
        "sdd"_a=sdd_energies,
        "bsa"_a=bsa_energies
    );

    results["interactions"] = py::dict(
        "atoms"_a=interaction_atoms,
        "sasa_atoms"_a=interaction_sasa_atoms,
        "contact_areas"_a=contact_areas,
        "distances"_a=distances,
        "dd_distances"_a=dd_distances,
        "dd_interacting"_a=dd_interacting
    );

    results["buried_area"] = py::dict(
        "atoms"_a=buried_by_atoms,
        "valid_atoms"_a=buried_by_atoms_valid,
        "areas"_a=buried_areas
    );

    results["overlaps"] = py::dict(
        "tables"_a=overlap_tables,
        "areas"_a=overlap_areas,
        "normalized_areas"_a=normalized_overlap_areas
    );

    if(include_matrix) {
        // Matrix handling remains the same
        map<vector<string>, vector<float>> matrixIJatom;
        map<vector<string>, vector<float>> matrixIJres;
        map<vector<string>, vector<string>> COLatom;
        map<vector<string>, vector<string>> COLres;
        map<vector<string>, vector<string>> ROWatom;
        map<vector<string>, vector<string>> ROWres;
        map<vector<string>, vector<uint32>> COLatomtype;
        map<vector<string>, vector<uint32>> ROWatomtype;

        GenerateInterBSAMatrix(atoms, matrixIJatom, matrixIJres,
                             COLatom, COLres, ROWatom, ROWres,
                             COLatomtype, ROWatomtype);

        results["matrices"] = py::dict(
            "atom_matrix"_a=matrixIJatom,
            "residue_matrix"_a=matrixIJres,
            "col_atoms"_a=COLatom,
            "row_atoms"_a=ROWatom,
            "col_types"_a=COLatomtype,
            "row_types"_a=ROWatomtype
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
    //PrintSASAResults(pdb,output);
    //PrintSplitAsaAtom(pdb,splitasa,atmasa_sasa);
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
    std::stringstream buffer;
    auto old_buf = std::cerr.rdbuf(buffer.rdbuf());
    auto atoms = PDBparser(pdb_file, "", true);
    if (atoms.empty()) {
        throw std::runtime_error("No atoms loaded from PDB file");
    }
    vdw_radii_.SetRadius(atoms, probe_radius_);
    std::cerr.rdbuf(old_buf);
    
    int Imode;
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
            cout << "#Protein-protein interface solver selected.\n";
            // For protein-protein, create chain groups from detected chains
            vector<string> chain_vec(proteinChains.begin(), proteinChains.end());
            chains = {vector<string>{chain_vec[0]}, vector<string>{chain_vec[1]}};
        } else {
            Imode = 4;  // Auto mode
            cout << "#Automatic interaction solver selected.\n";
        }
    } else {
        Imode = 1;  // Manual mode
        cout << "#Manual chain interaction solver selected.\n";
    }

    ChainSelector(chains, atoms);
    Generic_Solver(atoms, vdw_radii_.Points, chains, Imode, cl_mode_);
    cout << "After Generic_Solver - atoms with interactions: " 
        << count_if(atoms.begin(), atoms.end(), 
            [](const atom_struct& a) { return !a.INTERACTION_P.empty(); }) << "\n";
    GeneratePairInteractionData(atoms);
    cout << "After GeneratePairInteractionData - atoms with contact areas: "
     << count_if(atoms.begin(), atoms.end(),
        [](const atom_struct& a) { return !a.CONTACT_AREA.empty(); }) << "\n";

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
             py::arg("include_matrix") = true, // remove this
             "Calculate simple SASA using mode 0.\n"
             "Args:\n"
             "    pdb_file: Path to PDB file\n"
             "    include_matrix: Whether to include interaction matrices (default: True)\n"
             "Returns dict with SASA values, basic statistics, and optional matrices.");

    // Generic dSASA
    py::class_<GenericSASA>(m, "GenericSASA")
        .def(py::init<float, int>(),
            py::arg("probe_radius") = 1.4f,
            py::arg("compute_mode") = 0,
            R"(Initialize GenericSASA calculator.
                
                Args:
                    probe_radius: Water probe radius in Angstroms (default: 1.4)
                    compute_mode: Computation mode (default: 0))")
        .def("calculate", 
            [](GenericSASA& self, 
                const string& pdb_file,
                const vector<vector<string>>& chains_in,
                bool include_matrix) {
                vector<vector<string>> chains = chains_in;  // Create mutable copy
                return self.calculate(pdb_file, chains, include_matrix);
            },
            py::arg("pdb_file"),
            py::arg("chains") = vector<vector<string>>(),
            py::arg("include_matrix") = true,
            R"(Calculate delta SASA between chain groups.

                Args:
                    pdb_file: Path to PDB file
                    chains: List of chain groups (e.g. [['A'], ['B']] for chains A vs B)
                        Empty list triggers automatic mode
                    include_matrix: Whether to include interaction matrices (default: True)
                
                Mode Selection:
                    - Mode 1 (Manual): When specific chains are provided
                    - Mode 4 (Auto): When chains list is empty
                    - Mode 5 (Protein-protein): Automatically selected for exactly 2 protein chains
                
                Returns:
                    dict containing:
                    - Surface area values (SASA, buried surface area)
                    - Interface information and interactions
                    - Contact matrices (if include_matrix=True)
                
                Example:
                    >>> calc = GenericSASA()
                    >>> # Automatic mode
                    >>> results = calc.calculate('protein.pdb')
                    >>> # Manual chain selection
                    >>> results = calc.calculate('complex.pdb', chains=[['A'], ['B']])
                )");
    
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