#include "utils.hpp"
#include "constants.hpp"

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
    py::dict atom_data;
    py::list residue_data;
    py::dict residue_index;
    
    std::map<std::tuple<std::string, std::string, int>, py::dict> residue_temp_data;
    
    for (size_t i = 0; i < atoms.size(); ++i) {
        const auto& atom = atoms[i];
        if (!atom.ACTIVE) continue;

        float base_asa = atom.SASA;
        auto components = classify_atom_asa(atom, base_asa);
        
        // Create atom data with all fields
        py::dict atom_info = py::dict(
            "name"_a=atom.NAME,
            "resname"_a=atom.RESN,
            "chain"_a=atom.CHAIN,
            "resid"_a=atom.RESI,
            "struct_type"_a=atom.STRUCT_TYPE,
            "coords"_a=py::make_tuple(atom.COORDS[0], atom.COORDS[1], atom.COORDS[2]),
            "sphere_area"_a=atom.AREA,
            "sasa"_a=atom.SASA,
            "dsasa"_a=atom.EXT1,
            "vdw"_a=atom.VDW,
            "polar"_a=atom.POLAR,
            "charge"_a=atom.CHARGE,
            
            // All ASA components from our helper
            "total_asa"_a=base_asa,
            "bb_asa"_a=components.bb_asa,
            "sc_asa"_a=components.sc_asa,
            "majorgroove_asa"_a=components.majorgroove_asa,
            "minorgroove_asa"_a=components.minorgroove_asa,
            "nogroove_asa"_a=components.nogroove_asa,
            "polar_asa"_a=components.polar_asa,
            "polar_bb_asa"_a=components.polar_bb_asa,
            "polar_sc_asa"_a=components.polar_sc_asa,
            "polar_majorgroove_asa"_a=components.polar_majorgroove_asa,
            "polar_minorgroove_asa"_a=components.polar_minorgroove_asa,
            "polar_nogroove_asa"_a=components.polar_nogroove_asa,
            "hyd_asa"_a=components.hyd_asa,
            "hyd_bb_asa"_a=components.hyd_bb_asa,
            "hyd_sc_asa"_a=components.hyd_sc_asa,
            "hyd_majorgroove_asa"_a=components.hyd_majorgroove_asa,
            "hyd_minorgroove_asa"_a=components.hyd_minorgroove_asa,
            "hyd_nogroove_asa"_a=components.hyd_nogroove_asa,
            "lig_asa"_a=components.lig_asa,
            "lig_polar_asa"_a=components.lig_polar_asa,
            "lig_hyd_asa"_a=components.lig_hyd_asa
        );

        // Add to atom_data
        atom_data[py::str(std::to_string(atom.ID))] = std::move(atom_info);

        // Aggregate residue data with contacts and overlaps
        auto res_key = std::make_tuple(atom.CHAIN, atom.RESN, atom.RESI);
        if (residue_temp_data.find(res_key) == residue_temp_data.end()) {
            auto it = STANDARD_SASA_VALUES.find(atom.RESN);
            float standard_sasa = it != STANDARD_SASA_VALUES.end() ? it->second : 0.0f;
            
            residue_temp_data[res_key] = py::dict(
                "chain"_a=atom.CHAIN,
                "resname"_a=atom.RESN,
                "resid"_a=atom.RESI,
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
            res["total_sasa"] = res["total_sasa"].cast<float>() + atom.SASA;
            res["total_area"] = res["total_area"].cast<float>() + atom.AREA;
            res["n_atoms"] = res["n_atoms"].cast<int>() + 1;
            
            // Update center
            auto center = res["center"].cast<py::tuple>();
            res["center"] = py::make_tuple(
                center[0].cast<float>() + atom.COORDS[0],
                center[1].cast<float>() + atom.COORDS[1],
                center[2].cast<float>() + atom.COORDS[2]
            );
        }

        // Add contacts to residue data
        auto& res = residue_temp_data[res_key];
        auto contacts = res["contacts"].cast<py::dict>();
        for (uint32_t pos : atom.INTERACTION_SASA_P) {
            const auto& other_atom = atoms[pos];
            contacts[py::str(std::to_string(other_atom.ID))] = py::dict(
                "struct_type"_a=other_atom.STRUCT_TYPE,
                "contact_area"_a=atom.CONTACT_AREA.at(other_atom.ID),
                "distance"_a=atom.DISTANCES.at(pos)
            );
        }
        res["contacts"] = contacts;

        // Add overlaps to residue data
        auto overlaps = res["overlaps"].cast<py::list>();
        for (size_t j = 0; j < atom.ov_table.size(); ++j) {
            const auto& overlap_set = atom.ov_table[j];
            std::vector<uint32_t> overlap_ids;
            for (uint32_t pos : overlap_set) {
                overlap_ids.push_back(atoms[pos].ID);
            }
            
            overlaps.append(py::dict(
                "atoms"_a=overlap_ids,
                "overlap_area"_a=atom.ov_table_area[j],
                "normalized_area"_a=atom.ov_norm_area[j],
                "buried_area"_a=atom.AREA_BURIED_BY_ATOM_area[j]
            ));
        }
        res["overlaps"] = overlaps;
    }

    // Finalize residue data and build index
    size_t index = 0;
    for (auto& [key, res] : residue_temp_data) {
        const auto& [chain, resname, resid] = key;
        int n_atoms = res["n_atoms"].cast<int>();
        float total_sasa = res["total_sasa"].cast<float>();
        float standard_sasa = res["standard_sasa"].cast<float>();
        
        // Calculate dSASA properly
        float dsasa = std::max(0.0f, standard_sasa - total_sasa);
        
        // Update center
        auto center = res["center"].cast<py::tuple>();
        res["center"] = py::make_tuple(
            center[0].cast<float>() / n_atoms,
            center[1].cast<float>() / n_atoms,
            center[2].cast<float>() / n_atoms
        );
        
        // Set final dSASA
        res["dsasa"] = dsasa;
        
        // Add to residue_data list
        residue_data.append(res);
        
        // Create unique residue identifier and add to index
        std::string residue_id = chain + "_" + resname + "_" + std::to_string(resid);
        residue_index[py::str(residue_id)] = index++;
    }

    // Create final structure with both list and index
    results["atom_data"] = atom_data;
    results["residue_data"] = residue_data;
    results["residue_index"] = residue_index;
    
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
                atom_i.SASA = atom_i.AREA - atom_i.EXT0;  // Calculate SASA from buried area
                atom_i.EXT1 = atom_i.EXT0;  // In decoupled mode, dSASA equals buried area
            }
        }

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