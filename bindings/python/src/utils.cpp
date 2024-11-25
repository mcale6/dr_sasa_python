#include "utils.hpp"
#include "constants.hpp"

/*
The map index error usually occurs when we try to access values in maps (like CONTACT_AREA or DISTANCES) without checking if the key exists
Generic_Solver:

Sets INTERACTION_SASA_P
Populates CONTACT_AREA for each interaction
Populates DISTANCES for each interaction

DecoupledSolver:

Sets AREA_BURIED_BY_ATOM_vector/area
Doesn't automatically populate CONTACT_AREA and DISTANCES
 */

ASAComponents classify_atom_asa(const atom_struct& atom, float base_asa) {
    ASAComponents components;
    
    if (atom.MOL_TYPE == "PROTEIN") {
        if (is_backbone(atom.NAME)) {
            components.bb_asa = base_asa;
            if (is_polar_backbone(atom.NAME)) {
                components.polar_asa = base_asa;
                components.polar_bb_asa = base_asa;
            } else if (is_hyd_backbone(atom.NAME)) {
                components.hyd_asa = base_asa;
                components.hyd_bb_asa = base_asa;
            }
        } else {
            components.sc_asa = base_asa;
            if (is_polar_sidechain(atom.NAME)) {
                components.polar_asa = base_asa;
                components.polar_sc_asa = base_asa;
            } else if (is_hyd_sidechain(atom.NAME)) {
                components.hyd_asa = base_asa;
                components.hyd_sc_asa = base_asa;
            }
        }
    } else if (atom.MOL_TYPE == "DNA" || atom.MOL_TYPE == "RNA") {
        if (is_backbone(atom.NAME)) {
            components.bb_asa = base_asa;
            if (is_polar_backbone(atom.NAME)) {
                components.polar_asa = base_asa;
                components.polar_bb_asa = base_asa;
            } else if (is_hyd_backbone(atom.NAME)) {
                components.hyd_asa = base_asa;
                components.hyd_bb_asa = base_asa;
            }
        } else if (is_in_major_groove(atom.RESN, atom.NAME)) {
            components.majorgroove_asa = base_asa;
            if (atom.POLAR) {
                components.polar_asa = base_asa;
                components.polar_majorgroove_asa = base_asa;
            } else {
                components.hyd_asa = base_asa;
                components.hyd_majorgroove_asa = base_asa;
            }
        } else if (is_in_minor_groove(atom.RESN, atom.NAME)) {
            components.minorgroove_asa = base_asa;
            if (atom.POLAR) {
                components.polar_asa = base_asa;
                components.polar_minorgroove_asa = base_asa;
            } else {
                components.hyd_asa = base_asa;
                components.hyd_minorgroove_asa = base_asa;
            }
        } else {
            components.nogroove_asa = base_asa;
            if (atom.POLAR) {
                components.polar_asa = base_asa;
                components.polar_nogroove_asa = base_asa;
            } else {
                components.hyd_asa = base_asa;
                components.hyd_nogroove_asa = base_asa;
            }
        }
    } else if (atom.MOL_TYPE == "LIGAND") {
        components.lig_asa = base_asa;
        if (atom.POLAR) {
            components.lig_polar_asa = base_asa;
        } else {
            components.lig_hyd_asa = base_asa;
        }
    }
    
    return components;
}

py::dict create_analysis_results(const std::vector<atom_struct>& atoms, bool include_matrix) {
    py::dict results;
    py::dict atoms_dict;
    py::dict residues_dict;
    py::dict chains_dict;
    py::dict lookup_dict;
    py::dict metadata;
    
    // Track chain data
    std::map<std::string, std::set<std::string>> chain_to_residues;
    std::map<std::string, std::string> chain_types;
    std::set<std::string> all_types;
    std::map<std::tuple<std::string, std::string, int>, py::dict> residue_temp_data;

    for (size_t i = 0; i < atoms.size(); ++i) {
        const auto& atom = atoms[i];
        if (!atom.ACTIVE) continue;
        float base_asa = atom.SASA;
        auto components = classify_atom_asa(atom, base_asa);

        // Create backbone components
        py::dict backbone_components = py::dict(
            "total"_a=components.bb_asa,
            "polar"_a=components.polar_bb_asa,
            "hydrophobic"_a=components.hyd_bb_asa
        );

        // Create sidechain components
        py::dict sidechain_components = py::dict(
            "total"_a=components.sc_asa,
            "polar"_a=components.polar_sc_asa,
            "hydrophobic"_a=components.hyd_sc_asa
        );

        // Create groove components
        py::dict groove_components = py::dict(
            "major"_a=py::dict(
                "total"_a=components.majorgroove_asa,
                "polar"_a=components.polar_majorgroove_asa,
                "hydrophobic"_a=components.hyd_majorgroove_asa
            ),
            "minor"_a=py::dict(
                "total"_a=components.minorgroove_asa,
                "polar"_a=components.polar_minorgroove_asa,
                "hydrophobic"_a=components.hyd_minorgroove_asa
            ),
            "none"_a=py::dict(
                "total"_a=components.nogroove_asa,
                "polar"_a=components.polar_nogroove_asa,
                "hydrophobic"_a=components.hyd_nogroove_asa
            )
        );

        // Create ligand components
        py::dict ligand_components = py::dict(
            "total"_a=components.lig_asa,
            "polar"_a=components.lig_polar_asa,
            "hydrophobic"_a=components.lig_hyd_asa
        );

        // Create hierarchical atom data
        py::dict surface_info = py::dict(
            // Basic metrics
            "sphere_area"_a=atom.AREA,
            "sasa"_a=atom.SASA,
            "buried_area"_a=atom.EXT0,
            "contact_area"_a=atom.EXT1,
            "dsasa"_a=atom.EXT1,
            "total_asa"_a=base_asa,
            
            // Organized components
            "backbone"_a=backbone_components,
            "sidechain"_a=sidechain_components,
            "groove"_a=groove_components,
            "ligand"_a=ligand_components
        );

        py::dict properties = py::dict(
            "vdw"_a=atom.VDW,
            "polar"_a=atom.POLAR,
            "charge"_a=atom.CHARGE,
            "struct_type"_a=atom.STRUCT_TYPE
        );

        // Create contacts section with safety checks
        py::dict nonbonded;
        for (const auto& [contact_id, area] : atom.CONTACT_AREA) {
            auto distance_iter = atom.DISTANCES.find(contact_id);
            if (distance_iter != atom.DISTANCES.end()) {
                nonbonded[py::str(std::to_string(contact_id))] = py::dict(
                    "area"_a=area,
                    "distance"_a=distance_iter->second
                );
            }
        }

        py::list overlap_groups;
        for (size_t j = 0; j < atom.ov_table.size(); ++j) {
            if (j < atom.ov_table_area.size() && j < atom.ov_norm_area.size() && 
                j < atom.AREA_BURIED_BY_ATOM_area.size()) {
                py::list overlap_ids;
                for (uint32_t id : atom.ov_table[j]) {
                    overlap_ids.append(py::str(std::to_string(id)));
                }
                overlap_groups.append(py::dict(
                    "atoms"_a=overlap_ids,
                    "area"_a=atom.ov_table_area[j],
                    "normalized_area"_a=atom.ov_norm_area[j],
                    "buried_area"_a=atom.AREA_BURIED_BY_ATOM_area[j]
                ));
            }
        }

        py::dict contacts = py::dict(
            "nonbonded"_a=nonbonded,
            "overlap_groups"_a=overlap_groups
        );

        // Combine all atom info
        atoms_dict[py::str(std::to_string(atom.ID))] = py::dict(
            "name"_a=atom.NAME,
            "resid"_a=atom.RESI,
            "resname"_a=atom.RESN,
            "chain"_a=atom.CHAIN,
            "index"_a=i,
            "coords"_a=py::make_tuple(atom.COORDS[0], atom.COORDS[1], atom.COORDS[2]),
            "surface"_a=surface_info,
            "properties"_a=properties,
            "contacts"_a=contacts
        );

        // Aggregate residue data
        auto res_key = std::make_tuple(atom.CHAIN, atom.RESN, atom.RESI);
        if (residue_temp_data.find(res_key) == residue_temp_data.end()) {
            auto it = STANDARD_SASA_VALUES.find(atom.RESN);
            float standard_sasa = it != STANDARD_SASA_VALUES.end() ? it->second : 0.0f;
            
            residue_temp_data[res_key] = py::dict(
                "chain"_a=atom.CHAIN,
                "name"_a=atom.RESN,
                "number"_a=atom.RESI,
                "total_sasa"_a=atom.SASA,
                "total_area"_a=atom.AREA,
                "standard_sasa"_a=standard_sasa,
                "n_atoms"_a=1,
                "center"_a=py::make_tuple(atom.COORDS[0], atom.COORDS[1], atom.COORDS[2]),
                "contacts"_a=py::dict(),
                "overlaps"_a=py::list()
            );
        } else {
            auto& res = residue_temp_data[res_key];
            float total_sasa = res["total_sasa"].cast<float>() + atom.SASA;
            float total_area = res["total_area"].cast<float>() + atom.AREA;
            int n_atoms = res["n_atoms"].cast<int>() + 1;
            auto center = res["center"].cast<py::tuple>();
            
            res["total_sasa"] = total_sasa;
            res["total_area"] = total_area;
            res["n_atoms"] = n_atoms;
            res["center"] = py::make_tuple(
                center[0].cast<float>() + atom.COORDS[0],
                center[1].cast<float>() + atom.COORDS[1],
                center[2].cast<float>() + atom.COORDS[2]
            );

            // Add contacts to residue data with safety checks
            auto contacts = res["contacts"].cast<py::dict>();
            for (uint32_t pos : atom.INTERACTION_SASA_P) {
                if (pos < atoms.size()) {  // Bounds check
                    const auto& other_atom = atoms[pos];
                    auto contact_iter = atom.CONTACT_AREA.find(other_atom.ID);
                    auto distance_iter = atom.DISTANCES.find(pos);
                    if (contact_iter != atom.CONTACT_AREA.end() && distance_iter != atom.DISTANCES.end()) {
                        contacts[py::str(std::to_string(other_atom.ID))] = py::dict(
                            "struct_type"_a=other_atom.STRUCT_TYPE,
                            "contact_area"_a=contact_iter->second,
                            "distance"_a=distance_iter->second
                        );
                    }
                }
            }
            res["contacts"] = contacts;

            // Add overlaps to residue data with safety checks
            auto overlaps = res["overlaps"].cast<py::list>();
            for (size_t j = 0; j < atom.ov_table.size(); ++j) {
                if (j < atom.ov_table_area.size() && j < atom.ov_norm_area.size() && 
                    j < atom.AREA_BURIED_BY_ATOM_area.size()) {
                    std::vector<uint32_t> overlap_ids;
                    for (uint32_t pos : atom.ov_table[j]) {
                        if (pos < atoms.size()) {  // Bounds check
                            overlap_ids.push_back(atoms[pos].ID);
                        }
                    }
                    if (!overlap_ids.empty()) {
                        overlaps.append(py::dict(
                            "atoms"_a=overlap_ids,
                            "overlap_area"_a=atom.ov_table_area[j],
                            "normalized_area"_a=atom.ov_norm_area[j],
                            "buried_area"_a=atom.AREA_BURIED_BY_ATOM_area[j]
                        ));
                    }
                }
            }
            res["overlaps"] = overlaps;
        }

        // Track chain information
        std::string res_id = atom.CHAIN + "_" + atom.RESN + "_" + std::to_string(atom.RESI);
        chain_to_residues[atom.CHAIN].insert(res_id);
        chain_types[atom.CHAIN] = atom.MOL_TYPE;
        all_types.insert(atom.MOL_TYPE);
    }

    // Process residues
    for (const auto& [key, res] : residue_temp_data) {
        const auto& [chain, resname, resid] = key;
        std::string res_id = chain + "_" + resname + "_" + std::to_string(resid);

        float total_sasa = res["total_sasa"].cast<float>();
        float total_area = res["total_area"].cast<float>();
        float standard_sasa = res["standard_sasa"].cast<float>();
        int n_atoms = res["n_atoms"].cast<int>();
        auto center = res["center"].cast<py::tuple>();
        
        // Calculate center of mass
        auto final_center = py::make_tuple(
            center[0].cast<float>() / n_atoms,
            center[1].cast<float>() / n_atoms,
            center[2].cast<float>() / n_atoms
        );
        
        // Calculate dSASA
        float dsasa = std::max(0.0f, standard_sasa - total_sasa);

        residues_dict[py::str(res_id)] = py::dict(
            "identifiers"_a=py::dict(
                "chain"_a=chain,
                "name"_a=resname,
                "number"_a=resid
            ),
            "structure"_a=py::dict(
                "center"_a=final_center,
                "n_atoms"_a=n_atoms
            ),
            "surface"_a=py::dict(
                "total_sasa"_a=total_sasa,
                "total_area"_a=total_area,
                "standard_sasa"_a=standard_sasa,
                "dsasa"_a=dsasa
            ),
            "contacts"_a=res["contacts"],
            "overlaps"_a=res["overlaps"]
        );
    }

    // Create chains data with py::str keys
    for (const auto& [chain, residues] : chain_to_residues) {
        py::list residue_list;
        for (const auto& res_id : residues) {
            residue_list.append(res_id);
        }
        chains_dict[py::str(chain)] = py::dict(
            "type"_a=chain_types[chain],
            "residues"_a=residue_list
        );
    }

    // Create lookup tables with py::str keys
    py::dict atom_lookup;
    py::dict residue_lookup;
    
    py::dict residue_by_chain;
    for (const auto& [chain, residues] : chain_to_residues) {
        py::list residue_list;
        for (const auto& res_id : residues) {
            residue_list.append(res_id);
        }
        residue_by_chain[py::str(chain)] = residue_list;
    }
    residue_lookup[py::str("by_chain")] = residue_by_chain;
    
    lookup_dict[py::str("atoms")] = atom_lookup;
    lookup_dict[py::str("residues")] = residue_lookup;

    // Add metadata
    py::list chain_list;
    py::list type_list;
    for (const auto& [chain, _] : chain_types) {
        chain_list.append(chain);
    }
    for (const auto& type : all_types) {
        type_list.append(type);
    }
    
    metadata["total_atoms"] = static_cast<int>(atoms.size());
    metadata["total_residues"] = static_cast<int>(residue_temp_data.size());
    metadata["chains"] = chain_list;
    metadata["types"] = type_list;
    metadata["probe_radius"] = 1.4f;

    // Combine everything
    results["atoms"] = atoms_dict;
    results["residues"] = residues_dict;
    results["chains"] = chains_dict;
    results["lookup"] = lookup_dict;
    results["metadata"] = metadata;

    return results;
}

py::dict generate_inter_bsa_matrices(std::vector<atom_struct>& atoms) {
    map<vector<string>, vector<float>> matrixIJatom;
    map<vector<string>, vector<float>> matrixIJres;
    map<vector<string>, vector<string>> COLatom, COLres;
    map<vector<string>, vector<string>> ROWatom, ROWres;
    map<vector<string>, vector<uint32_t>> COLatomtype, ROWatomtype;

    GenerateInterBSAMatrix(atoms, matrixIJatom, matrixIJres, 
                          COLatom, COLres, ROWatom, ROWres,
                          COLatomtype, ROWatomtype);

    py::dict results;
    
    // Add atom-level matrices
    for (const auto& [key, matrix] : matrixIJatom) {
        const auto& cols = COLatom[key];
        const auto& rows = ROWatom[key];
        
        results[py::str(key[0] + "_vs_" + key[1])] = py::dict(
            "matrix"_a = py::array_t<float>({rows.size(), cols.size()}, matrix.data()),
            "row_labels"_a = rows,
            "col_labels"_a = cols,
            "row_types"_a = ROWatomtype[key],
            "col_types"_a = COLatomtype[key]
        );
    }
    
    // Add residue-level matrices
    py::dict res_matrices;
    for (const auto& [key, matrix] : matrixIJres) {
        const auto& cols = COLres[key];
        const auto& rows = ROWres[key];
        
        res_matrices[py::str(key[0] + "_vs_" + key[1])] = py::dict(
            "matrix"_a = py::array_t<float>({rows.size(), cols.size()}, matrix.data()),
            "row_labels"_a = rows,
            "col_labels"_a = cols
        );
    }
    results["residue_matrices"] = res_matrices;
    
    return results;
}

py::dict generate_intra_bsa_matrices(std::vector<atom_struct>& atoms) {
    vector<float> matrixIJatom, matrixIJres;
    vector<string> COLatom, COLres, ROWatom, ROWres;
    vector<uint32_t> COLatomtype, ROWatomtype;

    GenerateIntraBSAMatrix(atoms, matrixIJatom, matrixIJres,
                          COLatom, COLres, ROWatom, ROWres,
                          COLatomtype, ROWatomtype);

    return py::dict(
        "atom_matrix"_a = py::array_t<float>(
            {ROWatom.size(), COLatom.size()}, 
            matrixIJatom.data()
        ),
        "residue_matrix"_a = py::array_t<float>(
            {ROWres.size(), COLres.size()}, 
            matrixIJres.data()
        ),
        "atom_labels"_a = COLatom,
        "residue_labels"_a = COLres,
        "atom_types"_a = COLatomtype
    );
}

void calculate_contact_areas_from_overlaps(vector<atom_struct>& pdb) {
    // Process each atom
    for (size_t pos = 0; pos < pdb.size(); ++pos) {
        auto& atom_i = pdb[pos];
        if (!atom_i.ACTIVE) continue;

        vector<uint32> o_atoms;
        vector<uint32> interac_pos;
        vector<char> valid_overlaps(atom_i.AREA_BURIED_BY_ATOM_vector.size(), true);

        // Find interacting atoms with different STRUCT_TYPE
        for (uint32 i = 0; i < atom_i.INTERACTION_P.size(); ++i) {
            if (atom_i.STRUCT_TYPE != pdb[atom_i.INTERACTION_P[i]].STRUCT_TYPE) {
                interac_pos.push_back(atom_i.INTERACTION_P[i]);
            }
        }

        // Validate overlaps and collect unique interacting atoms
        for (uint32 i = 0; i < atom_i.AREA_BURIED_BY_ATOM_vector.size(); ++i) {
            for (uint32 k = 0; k < atom_i.AREA_BURIED_BY_ATOM_vector[i].size(); ++k) {
                uint32 other = atom_i.AREA_BURIED_BY_ATOM_vector[i][k];
                if (interac_pos.end() == find(interac_pos.begin(), interac_pos.end(), other)) {
                    valid_overlaps[i] = false;
                }
            }
            if (!valid_overlaps[i]) continue;

            atom_i.AREA_BURIED_BY_ATOM_vector_valid.push_back(i);
            for (uint32 j = 0; j < atom_i.AREA_BURIED_BY_ATOM_vector[i].size(); ++j) {
                if (o_atoms.end() == find(o_atoms.begin(), o_atoms.end(), atom_i.AREA_BURIED_BY_ATOM_vector[i][j])) {
                    o_atoms.push_back(atom_i.AREA_BURIED_BY_ATOM_vector[i][j]);
                }
            }
        }

        sort(o_atoms.begin(), o_atoms.end());

        // Calculate EXT0 (total buried area)
        atom_i.EXT0 = 0.0;
        for (uint32 i = 0; i < atom_i.AREA_BURIED_BY_ATOM_area.size(); ++i) {
            if (valid_overlaps[i]) {
                atom_i.EXT0 += atom_i.AREA_BURIED_BY_ATOM_area[i];
            }
        }
        // because decoupled
        atom_i.SASA = atom_i.AREA - atom_i.EXT0;  // SASA is what's not buried
        atom_i.EXT1 = atom_i.EXT0;                // dSASA equals buried area in decoupled

        // Store interaction atoms and calculate contact areas
        atom_i.INTERACTION_SASA_P = o_atoms;
        
        // Calculate contact areas
        for (uint32 other_pos : o_atoms) {
            float T_area = 0;
            auto& atom_j = pdb[other_pos];
            
            for (uint32 j = 0; j < atom_i.AREA_BURIED_BY_ATOM_vector.size(); ++j) {
                if (valid_overlaps[j] && 
                    atom_i.AREA_BURIED_BY_ATOM_vector[j].end() != 
                    find(atom_i.AREA_BURIED_BY_ATOM_vector[j].begin(), 
                         atom_i.AREA_BURIED_BY_ATOM_vector[j].end(), 
                         other_pos)) {
                    T_area += atom_i.AREA_BURIED_BY_ATOM_area[j];
                }
            }
            
            atom_i.CONTACT_AREA[atom_j.ID] = T_area;
        }

        // Populate overlap tables
        atom_i.ov_table.clear();
        atom_i.ov_table_area.clear();
        atom_i.ov_norm_area.clear();

        for (uint32 i = 0; i < atom_i.AREA_BURIED_BY_ATOM_vector.size(); ++i) {
            if (valid_overlaps[i]) {
                atom_i.ov_table.push_back(atom_i.AREA_BURIED_BY_ATOM_vector[i]);
                atom_i.ov_table_area.push_back(atom_i.AREA_BURIED_BY_ATOM_area[i]);
                atom_i.ov_norm_area.push_back(atom_i.AREA_BURIED_BY_ATOM_area[i] / 
                                            atom_i.AREA_BURIED_BY_ATOM_vector[i].size());
            }
        }
    }
}