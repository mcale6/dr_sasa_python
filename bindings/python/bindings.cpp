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
//#include "interaction_analysis.h"  // Add the new header

namespace py = pybind11;
using namespace py::literals;
//using namespace py::drsasa;

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
extern void CalculateDNA_ProtInteractions(vector<atom_struct>& pdb, int mode);

// Additional useful external functions we might need later
// GenerateInterBSAMatrix 
// GenerateIntraBSAMatrix 
// GeneratePDistribution 

// differences? add per atom interactions
// PrintDNA_ProtResults for overlaps from CalculateDNA_ProtInteractions DNA /Prtoein?
// PrintDNA_ProtResultsByAtomMatrix matrix from CalculateDNA_ProtInteractions
// PrintDNA_ProtInteractions ?

// - SimpleSASA      (Mode 0)  # Basic SASA calculation
// - GenericSASA     (Mode 1)  # Chain/molecular type interactions
// - DecoupledSASA   (Mode 4)  # Decoupled interface analysis
// - RelativeSASA    (Mode 100) # Relative accessibility this is wrong in the SolverDataProcessing per atom not per residue!
 
// "Missing"
// - Mode 2 (Residue dSASA) # returned in the dict handled by GenericSASA/DecoupledSASA
// - Mode 3 (Atom dSASA)    # returned in the dict andled by GenericSASA/DecoupledSASA

static constexpr float DEFAULT_PROBE_RADIUS = 1.4f;  // Water probe in Angstroms

struct SurfaceSummary {
    float complex_surface;      // Total surface in complex (same for both)
    float buried_by_other;      // Surface buried by other chain
    float buries_in_other;      // Surface this chain buries in other
    float interface_area;       // Average of both buried surfaces
    string object_id;          // Chain/object identifier
};

// Helper function to convert matrix data to numpy array with metadata
py::dict create_matrix_dict(const vector<float>& matrix_data,
                          const vector<string>& col_labels,
                          const vector<string>& row_labels,
                          size_t num_rows,
                          size_t num_cols) {
    // Create numpy array with proper shape
    vector<ssize_t> shape = {(ssize_t)num_rows, (ssize_t)num_cols};
    vector<ssize_t> strides = {(ssize_t)(num_cols * sizeof(float)), sizeof(float)};
    
    return py::dict(
        "values"_a=py::array_t<float>(shape, strides, matrix_data.data()),
        "row_labels"_a=row_labels,
        "col_labels"_a=col_labels
    );
}

py::dict calculate_interaction_data(vector<atom_struct>& pdb) {
    // Get unique chains/objects
    set<string> object_types;
    for (auto& atom : pdb) {
        object_types.insert(atom.STRUCT_TYPE);
    }
    vector<string> objs(object_types.begin(), object_types.end());

    // Initialize surfaces tracking
    float total_complex_surface = 0.0f;
    map<pair<string, string>, float> buried_between_objects;
    
    // Build indices just like Print_MatrixInsideAtom
    vector<string> atom_labels;
    vector<string> atom_full_labels;
    map<string, uint32_t> atom_map;
    vector<string> res_labels;
    map<string, uint32_t> res_map;

    // First pass - collect labels and complex surface
    for (auto& atom : pdb) {
        total_complex_surface += atom.SASA;
        string aID = atom.sID();
        string rID = atom.rsID();
        
        stringstream atomtype;
        atomtype << aID << "/" << atom.MOL_TYPE << "/" << atom.ATOM_TYPE;
        
        atom_labels.push_back(aID);
        atom_full_labels.push_back(atomtype.str());
        atom_map[aID] = atom_labels.size() - 1;
        res_labels.push_back(rID);
    }

    // Process residue labels
    sort(res_labels.begin(), res_labels.end());
    res_labels.erase(unique(res_labels.begin(), res_labels.end()), res_labels.end());
    for (uint32_t k = 0; k < res_labels.size(); ++k) {
        res_map[res_labels[k]] = k;
    }

    // Initialize matrices -> PrintDNA_ProtResultsByAtomMatrix
    size_t num_atoms = atom_labels.size();
    size_t num_residues = res_labels.size();
    vector<float> atom_matrix(num_atoms * num_atoms, NAN);
    vector<float> res_matrix(num_residues * num_residues, NAN);

    for (uint32_t i = 0; i < pdb.size(); ++i) {
        auto& atom = pdb[i];
        string aID = atom.sID();
        string rID = atom.rsID();

        for (uint32_t pos : atom.INTERACTION_SASA_P) {
            if (atom.STRUCT_TYPE != pdb[pos].STRUCT_TYPE) {
                uint32_t col = atom_map[aID];
                uint32_t line = atom_map[pdb[pos].sID()];
                uint32_t col_res = res_map[rID];
                uint32_t line_res = res_map[pdb[pos].rsID()];

                uint64_t idx_atom = col + line * num_atoms;
                uint64_t idx_res = col_res + line_res * num_residues;

                if (std::isnan(atom_matrix[idx_atom])) atom_matrix[idx_atom] = 0.0;
                if (std::isnan(res_matrix[idx_res])) res_matrix[idx_res] = 0.0;

                for (uint32_t r = 0; r < atom.ov_table.size(); ++r) {
                    const auto& ov = atom.ov_table[r];
                    if (find(ov.begin(), ov.end(), pos) != ov.end()) {
                        float area = atom.ov_norm_area[r];
                        atom_matrix[idx_atom] += area;
                        res_matrix[idx_res] += area;
                        buried_between_objects[{atom.STRUCT_TYPE, pdb[pos].STRUCT_TYPE}] += area;
                    }
                }
            }
        }
    }
    // Add residue interaction analysis
    vector<pair<int,int>> residue_pairs;
    vector<float> buried_areas;
    vector<float> dsasa_res1;
    vector<float> dsasa_res2;
    
    // Track processed pairs and residue dSASA
    vector<pair<uint32,uint32>> processed_pairs;
    map<int, float> residue_dsasa;
    map<pair<int,int>, float> residue_contacts;

    // Calculate residue dSASA and contacts
    for (uint32_t i = 0; i < pdb.size(); ++i) {
        auto& atom_i = pdb[i];
        if (!atom_i.ACTIVE || atom_i.EXT1 == 0) continue;
        
        // Add to residue's total dSASA
        if (residue_dsasa.count(atom_i.RESI) == 0) {
            residue_dsasa[atom_i.RESI] = 0.0;
        }
        residue_dsasa[atom_i.RESI] += atom_i.EXT1;

        // Process interactions
        for (uint32_t j = 0; j < atom_i.INTERACTION_SASA_P.size(); ++j) {
            auto j_pos = atom_i.INTERACTION_SASA_P[j];
            auto& atom_j = pdb[j_pos];
            
            pair<uint32_t,uint32_t> p1(i, j_pos);
            pair<uint32_t,uint32_t> p2(j_pos, i);
            
            if (!atom_j.ACTIVE || atom_j.EXT1 == 0) continue;
            if (find(processed_pairs.begin(), processed_pairs.end(), p1) != processed_pairs.end() ||
                find(processed_pairs.begin(), processed_pairs.end(), p2) != processed_pairs.end()) {
                continue;
            }

            processed_pairs.push_back(p1);
            processed_pairs.push_back(p2);

            // Get contact areas in both directions
            float ji = atom_i.CONTACT_AREA.count(atom_j.ID) ? 
                      atom_i.CONTACT_AREA.at(atom_j.ID) : 0.0;
            float ij = atom_j.CONTACT_AREA.count(atom_i.ID) ? 
                      atom_j.CONTACT_AREA.at(atom_i.ID) : 0.0;

            // Store residue contact data
            auto res_ij = make_pair(atom_i.RESI, atom_j.RESI);
            auto res_ji = make_pair(atom_j.RESI, atom_i.RESI);
            
            if (residue_contacts.count(res_ij) == 0) residue_contacts[res_ij] = 0.0;
            if (residue_contacts.count(res_ji) == 0) residue_contacts[res_ji] = 0.0;
            
            residue_contacts[res_ij] += ji;
            residue_contacts[res_ji] += ij;
        }
    }

    // Convert to vectors for Python
    for (const auto& contact : residue_contacts) {
        const auto& res_pair = contact.first;
        float buried_area = contact.second;
        
        residue_pairs.push_back(res_pair);
        buried_areas.push_back(buried_area);
        dsasa_res1.push_back(residue_dsasa[res_pair.first]);
        dsasa_res2.push_back(residue_dsasa[res_pair.second]);
    }

    // Create residue interaction dictionary
    py::dict residue_interactions = py::dict(
        "residue1"_a=py::list(),
        "residue2"_a=py::list(),
        "buried_area"_a=py::list(),
        "dsasa_res1"_a=py::list(),
        "dsasa_res2"_a=py::list()
    );

    // Fill the lists
    for (size_t i = 0; i < residue_pairs.size(); ++i) {
        residue_interactions["residue1"].cast<py::list>().append(residue_pairs[i].first);
        residue_interactions["residue2"].cast<py::list>().append(residue_pairs[i].second);
        residue_interactions["buried_area"].cast<py::list>().append(buried_areas[i]);
        residue_interactions["dsasa_res1"].cast<py::list>().append(dsasa_res1[i]);
        residue_interactions["dsasa_res2"].cast<py::list>().append(dsasa_res2[i]);
    }

    // Atom interaction analysis
    vector<pair<int,int>> atom_pairs;
    vector<float> atom_buried_areas;
    vector<float> dsasa_atom1;
    vector<float> dsasa_atom2;

    // Track processed atom pairs
    vector<pair<uint32,uint32>> processed_atom_pairs;
    map<int, float> atom_dsasa;
    map<pair<int,int>, float> atom_contacts;

    // Calculate atom dSASA and contacts
    for (uint32_t i = 0; i < pdb.size(); ++i) {
        auto& atom_i = pdb[i];
        if (!atom_i.ACTIVE || atom_i.EXT1 == 0) continue;
        
        // Add to atom's total dSASA
        atom_dsasa[atom_i.ID] = atom_i.EXT1;

        // Process interactions
        for (uint32_t j = 0; j < atom_i.INTERACTION_SASA_P.size(); ++j) {
            auto j_pos = atom_i.INTERACTION_SASA_P[j];
            auto& atom_j = pdb[j_pos];
            
            pair<uint32_t,uint32_t> p1(i, j_pos);
            pair<uint32_t,uint32_t> p2(j_pos, i);
            
            if (!atom_j.ACTIVE || atom_j.EXT1 == 0) continue;
            if (find(processed_atom_pairs.begin(), processed_atom_pairs.end(), p1) != processed_atom_pairs.end() ||
                find(processed_atom_pairs.begin(), processed_atom_pairs.end(), p2) != processed_atom_pairs.end()) {
                continue;
            }

            processed_atom_pairs.push_back(p1);
            processed_atom_pairs.push_back(p2);

            // Get contact areas in both directions
            float ji = atom_i.CONTACT_AREA.count(atom_j.ID) ? 
                    atom_i.CONTACT_AREA.at(atom_j.ID) : 0.0;
            float ij = atom_j.CONTACT_AREA.count(atom_i.ID) ? 
                    atom_j.CONTACT_AREA.at(atom_i.ID) : 0.0;

            // Store atom contact data
            auto atom_pair_ij = make_pair(atom_i.ID, atom_j.ID);
            auto atom_pair_ji = make_pair(atom_j.ID, atom_i.ID);
            
            if (atom_contacts.count(atom_pair_ij) == 0) atom_contacts[atom_pair_ij] = 0.0;
            if (atom_contacts.count(atom_pair_ji) == 0) atom_contacts[atom_pair_ji] = 0.0;
            
            atom_contacts[atom_pair_ij] += ji;
            atom_contacts[atom_pair_ji] += ij;
        }
    }

    // Convert to vectors for Python
    for (const auto& contact : atom_contacts) {
        const auto& atom_pair = contact.first;
        float buried_area = contact.second;
        
        atom_pairs.push_back(atom_pair);
        atom_buried_areas.push_back(buried_area);
        dsasa_atom1.push_back(atom_dsasa[atom_pair.first]);
        dsasa_atom2.push_back(atom_dsasa[atom_pair.second]);
    }

    // Create atom interaction dictionary
    py::dict atom_interactions = py::dict(
        "atom1_id"_a=py::list(),
        "atom2_id"_a=py::list(),
        "buried_area"_a=py::list(),
        "dsasa_atom1"_a=py::list(),
        "dsasa_atom2"_a=py::list()
    );

    // Fill the lists
    for (size_t i = 0; i < atom_pairs.size(); ++i) {
        atom_interactions["atom1_id"].cast<py::list>().append(atom_pairs[i].first);
        atom_interactions["atom2_id"].cast<py::list>().append(atom_pairs[i].second);
        atom_interactions["buried_area"].cast<py::list>().append(atom_buried_areas[i]);
        atom_interactions["dsasa_atom1"].cast<py::list>().append(dsasa_atom1[i]);
        atom_interactions["dsasa_atom2"].cast<py::list>().append(dsasa_atom2[i]);
    }

    // Create summaries
    py::dict summary_dict;
    for (const auto& obj_A : objs) {
        for (const auto& obj_B : objs) {
            if (obj_A != obj_B) {
                float buried_by_other = buried_between_objects[{obj_A, obj_B}];
                float buries_in_other = buried_between_objects[{obj_B, obj_A}];
                float interface_area = (buried_by_other + buries_in_other) / 2.0f;

                summary_dict[py::str(obj_A)] = py::dict(
                    "complex_surface"_a=total_complex_surface,
                    "buried_by_other"_a=buried_by_other,     // A<---B
                    "buries_in_other"_a=buries_in_other,     // A--->B
                    "interface_area"_a=interface_area         // Average
                );
            }
        }
    }

    return py::dict(
        "matrices"_a=py::dict(
            "atom"_a=create_matrix_dict(atom_matrix, atom_full_labels, atom_full_labels, 
                                      num_atoms, num_atoms),
            "residue"_a=create_matrix_dict(res_matrix, res_labels, res_labels, 
                                         num_residues, num_residues)
        ),
        "surface_summary"_a=summary_dict,
        "residue_interactions"_a=residue_interactions,
        "atom_interactions"_a=atom_interactions

    );
}

py::dict create_analysis_results(vector<atom_struct>& atoms, bool include_matrix = false) {
    py::dict results;
    const size_t n_atoms = atoms.size();

    // Basic atom info arrays
    vector<string> residues, chains, elements;
    vector<string> mol_types;
    vector<int> residue_nums, atom_ids;
    vector<bool> is_hetatm;
    vector<int> atom_types;
    vector<int> polarity;

    // Geometric properties
    vector<float> coordinates, radii, radii2, vdw;
    vector<float> b_factors;
    vector<string> charges;
    
    // Surface areas and volumes
    auto sasa = py::array_t<double>(n_atoms);
    auto area = py::array_t<double>(n_atoms);
    auto ext0 = py::array_t<double>(n_atoms);
    auto ext1 = py::array_t<double>(n_atoms);
    auto accessible = py::array_t<double>(n_atoms);
    
    // Interaction data 
    vector<vector<uint32_t>> interaction_atoms;
    vector<vector<uint32_t>> interaction_sasa_atoms;
    vector<map<uint32_t, float>> contact_areas;
        
    // Buried area data
    vector<vector<vector<uint32_t>>> buried_by_atoms;
    vector<vector<uint32_t>> buried_by_atoms_valid;
    vector<vector<float>> buried_areas;
        
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
            atom_ids.push_back(i);
            residues.push_back(atom.RESN);
            chains.push_back(atom.CHAIN);
            elements.push_back(atom.ELEMENT);
            mol_types.push_back(atom.MOL_TYPE);
            residue_nums.push_back(atom.RESI);
            // Flags
            is_hetatm.push_back(atom.HETATM);
            // Types
            atom_types.push_back(atom.ATOM_TYPE);
            polarity.push_back(atom.POLAR);
            // Geometric properties
            coordinates.insert(coordinates.end(), atom.COORDS.begin(), atom.COORDS.end());
            radii.push_back(atom.RADIUS);
            radii2.push_back(atom.RADIUS2);
            vdw.push_back(atom.VDW);
            b_factors.push_back(atom.TFACTOR);
            // Surface areas
            sasa_ptr[i] = atom.SASA;
            area_ptr[i] = atom.AREA;
            ext0_ptr[i] = atom.EXT0;
            ext1_ptr[i] = atom.EXT1;
            acc_ptr[i] = atom.ACCESSIBLE_P;

            // Interaction data
            interaction_atoms.push_back(atom.INTERACTION_P);
            interaction_sasa_atoms.push_back(atom.INTERACTION_SASA_P);
            contact_areas.push_back(atom.CONTACT_AREA);

            // Buried area data
            buried_by_atoms.push_back(atom.AREA_BURIED_BY_ATOM_vector);
            buried_by_atoms_valid.push_back(atom.AREA_BURIED_BY_ATOM_vector_valid);
            buried_areas.push_back(atom.AREA_BURIED_BY_ATOM_area);
            
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
        "ids"_a=atom_ids,
        "residues"_a=residues,
        "chains"_a=chains,
        "elements"_a=elements,
        "mol_types"_a=mol_types,
        "residue_numbers"_a=residue_nums,
        "is_hetatm"_a=is_hetatm,
        "atom_types"_a=atom_types,
        "polarity"_a=polarity,
        "coordinates"_a=py::array_t<float>(coord_shape, coord_strides, coordinates.data()),
        "radii"_a=radii,
        "radii2"_a=radii2,
        "vdw"_a=vdw,
        "b_factors"_a=b_factors
    );

    results["surface"] = py::dict( // not clear needs fix
        "area"_a=area,
        "sasa"_a=sasa,
        "buried_area"_a=ext0,
        "fast_dsasa"_a=ext1,
        "accessible"_a=accessible
    );

    results["interactions"] = py::dict(
        "atoms"_a=interaction_atoms,
        "sasa_atoms"_a=interaction_sasa_atoms,
        "contact_areas"_a=contact_areas
    );

    results["buried_area"] = py::dict(
        "atoms"_a=buried_by_atoms,
        "valid_atoms"_a=buried_by_atoms_valid,
        "areas"_a=buried_areas
    );

    results["overlaps"] = py::dict( // excluded in uttils
        "tables"_a=overlap_tables,
        "areas"_a=overlap_areas,
        "normalized_areas"_a=normalized_overlap_areas
    );

    if (include_matrix) {
        // Calculate interactions data once
        auto interaction_data = calculate_interaction_data(atoms);
        
        // Extract both matrices and summary from the same calculation
        results["matrices"] = interaction_data["matrices"];
        results["surface_summary"] = interaction_data["surface_summary"];
        results["residue_interactions"] = interaction_data["residue_interactions"];
        results["atom_interactions"] = interaction_data["atom_interactions"];

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
    
    py::dict calculate(const string& pdb_file);
};

py::dict SimpleSASA::calculate(const string& pdb_file) {
    std::stringstream buffer;
    auto old_buf = std::cerr.rdbuf(buffer.rdbuf());
    auto atoms = PDBparser(pdb_file, "", true);
    if (atoms.empty()) {
        throw std::runtime_error("No atoms loaded from PDB file");
    }
    vdw_radii_.SetRadius(atoms, probe_radius_);
    std::cerr.rdbuf(old_buf);

    SolveInteractions(atoms, 0);  // Mode 0 specific
    SimpleSolverCL(atoms, vdw_radii_.Points, cl_mode_);
    //PrintSASAResults(pdb,output);
    //PrintSplitAsaAtom(pdb,splitasa,atmasa_sasa);
    return create_analysis_results(atoms); // include matrix false default
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
    //cout << "After Generic_Solver - atoms with interactions: " 
    //    << count_if(atoms.begin(), atoms.end(), 
    //        [](const atom_struct& a) { return !a.INTERACTION_P.empty(); }) << "\n";
    GeneratePairInteractionData(atoms);
    CalculateDNA_ProtInteractions(atoms, cl_mode_);
    //cout << "After GeneratePairInteractionData - atoms with contact areas: "
    // << count_if(atoms.begin(), atoms.end(),
    //    [](const atom_struct& a) { return !a.CONTACT_AREA.empty(); }) << "\n";
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

// RelSASA class is wrong and removed

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
             "Args:\n"
             "    pdb_file: Path to PDB file\n"
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

    // Convenience functions
    m.def("calculate_simple_sasa", 
        [](const string& pdb_file, float probe_radius = 1.4f) {
            SimpleSASA calculator(probe_radius);
            return calculator.calculate(pdb_file);
        }, 
        py::arg("pdb_file"), 
        py::arg("probe_radius") = 1.4f,
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