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
    py::dict atoms_dict;
    py::dict residues_dict;
    py::dict chains_dict;
    py::dict lookup;
    std::set<std::string> mol_types;

    // Initialize lookup tables
    py::dict atom_lookup;
    py::dict res_lookup;
    atom_lookup["by_id"] = py::dict();
    atom_lookup["by_residue"] = py::dict();
    atom_lookup["by_chain"] = py::dict();
    res_lookup["by_chain"] = py::dict();
    
    // Track unique chains and build lookup tables
    std::map<std::string, std::vector<std::string>> chain_to_atoms;
    std::map<std::string, std::vector<std::string>> chain_to_residues;
    std::map<std::string, std::vector<std::string>> residue_to_atoms;
    
    // First pass: Process atoms and build initial lookup tables
    for (size_t i = 0; i < atoms.size(); ++i) {
        const auto& atom = atoms[i];
        if (!atom.ACTIVE) continue;

        // Ensure chain ID is properly handled
        std::string chain = atom.CHAIN;
        if (chain.empty()) {
            chain = "_";
        }

        mol_types.insert(atom.MOL_TYPE);
        std::string atom_id = std::to_string(atom.ID);
        std::string residue_id = chain + "_" + atom.RESN + "_" + std::to_string(atom.RESI);
        
        // Process contacts
        py::dict nonbonded;
        for (uint32_t pos : atom.INTERACTION_SASA_P) {
            if (pos >= atoms.size()) continue;

            const auto& other_atom = atoms[pos];
            std::string other_id = std::to_string(other_atom.ID);
            
            if (atom.CONTACT_AREA.count(other_atom.ID) && atom.DISTANCES.count(pos)) {
                nonbonded[py::str(other_id)] = py::dict(
                    "area"_a=atom.CONTACT_AREA.at(other_atom.ID),
                    "distance"_a=atom.DISTANCES.at(pos)
                );
            }
        }
        
        // Process overlaps
        py::list overlap_groups;
        for (size_t j = 0; j < atom.ov_table.size(); ++j) {
            if (j >= atom.AREA_BURIED_BY_ATOM_area.size()) break;

            const auto& overlap_set = atom.ov_table[j];
            std::vector<std::string> overlap_atoms;
            
            bool valid_overlap = true;
            for (uint32_t pos : overlap_set) {
                if (pos >= atoms.size()) {
                    valid_overlap = false;
                    break;
                }
                overlap_atoms.push_back(std::to_string(atoms[pos].ID));
            }
            
            if (!valid_overlap) continue;
            
            float buried_area = atom.AREA_BURIED_BY_ATOM_area[j];
            float norm_area = overlap_set.empty() ? 0.0f : buried_area / overlap_set.size();

            if (j < atom.ov_table_area.size()) {
                overlap_groups.append(py::dict(
                    "atoms"_a=overlap_atoms,
                    "area"_a=atom.ov_table_area[j],
                    "normalized_area"_a=norm_area,
                    "buried_area"_a=buried_area
                ));
            }
        }

        // Calculate surface values
        float total_buried_area = 0.0f;
        for (const auto& area : atom.AREA_BURIED_BY_ATOM_area) {
            total_buried_area += area;
        }

        float total_contact_area = 0.0f;
        for (const auto& [_, area] : atom.CONTACT_AREA) {
            total_contact_area += area;
        }

        // Create atom entry
        atoms_dict[py::str(atom_id)] = py::dict(
            // Basic Properties
            "name"_a=atom.NAME,
            "resid"_a=atom.RESI,
            "resname"_a=atom.RESN,
            "chain"_a=chain,
            "index"_a=i,
            "coords"_a=py::make_tuple(atom.COORDS[0], atom.COORDS[1], atom.COORDS[2]),
            
            // Surface Analysis
            "surface"_a=py::dict(
                "sphere_area"_a=atom.AREA,
                "sasa"_a=atom.SASA,
                "buried_area"_a=total_buried_area,
                "contact_area"_a=total_contact_area,
                "dsasa"_a=atom.EXT1
            ),
            
            // Physical Properties
            "properties"_a=py::dict(
                "vdw"_a=atom.VDW,
                "polar"_a=atom.POLAR,
                "charge"_a=atom.CHARGE,
                "struct_type"_a=atom.MOL_TYPE
            ),
            
            // Contacts
            "contacts"_a=py::dict(
                "bonded"_a=py::list(),
                "nonbonded"_a=nonbonded,
                "overlap_groups"_a=overlap_groups
            )
        );
        
        // Update lookup tables
        chain_to_atoms[chain].push_back(atom_id);
        residue_to_atoms[residue_id].push_back(atom_id);
    }
    
    // Second pass: Build residue information
    for (const auto& [residue_id, atom_ids] : residue_to_atoms) {
        if (atom_ids.empty()) continue;

        float total_sasa = 0.0f;
        float total_area = 0.0f;
        std::vector<float> center = {0.0f, 0.0f, 0.0f};
        
        // Get first atom to extract residue info
        const auto& first_atom = atoms[std::stoi(atom_ids[0]) - 1];
        
        // Calculate standard SASA
        auto it = STANDARD_SASA_VALUES.find(first_atom.RESN);
        float standard_sasa = it != STANDARD_SASA_VALUES.end() ? it->second : 0.0f;
        
        // Initialize ASA components
        ASAComponents total_components;
        
        // Aggregate atom information
        for (const auto& atom_id : atom_ids) {
            const auto& atom = atoms[std::stoi(atom_id) - 1];
            total_sasa += atom.SASA;
            total_area += atom.AREA;
            
            center[0] += atom.COORDS[0];
            center[1] += atom.COORDS[1];
            center[2] += atom.COORDS[2];
            
            // Accumulate ASA components
            auto components = classify_atom_asa(atom, atom.SASA);
            total_components.bb_asa += components.bb_asa;
            total_components.sc_asa += components.sc_asa;
            total_components.majorgroove_asa += components.majorgroove_asa;
            total_components.minorgroove_asa += components.minorgroove_asa;
            total_components.nogroove_asa += components.nogroove_asa;
            total_components.polar_asa += components.polar_asa;
            total_components.polar_bb_asa += components.polar_bb_asa;
            total_components.polar_sc_asa += components.polar_sc_asa;
            total_components.polar_majorgroove_asa += components.polar_majorgroove_asa;
            total_components.polar_minorgroove_asa += components.polar_minorgroove_asa;
            total_components.polar_nogroove_asa += components.polar_nogroove_asa;
            total_components.hyd_asa += components.hyd_asa;
            total_components.hyd_bb_asa += components.hyd_bb_asa;
            total_components.hyd_sc_asa += components.hyd_sc_asa;
            total_components.hyd_majorgroove_asa += components.hyd_majorgroove_asa;
            total_components.hyd_minorgroove_asa += components.hyd_minorgroove_asa;
            total_components.hyd_nogroove_asa += components.hyd_nogroove_asa;
            total_components.lig_asa += components.lig_asa;
            total_components.lig_polar_asa += components.lig_polar_asa;
            total_components.lig_hyd_asa += components.lig_hyd_asa;
        }
        
        // Average the center
        float n_atoms = static_cast<float>(atom_ids.size());
        for (auto& coord : center) {
            coord /= n_atoms;
        }
        
        residues_dict[py::str(residue_id)] = py::dict(
            "identifiers"_a=py::dict(
                "chain"_a=first_atom.CHAIN,
                "name"_a=first_atom.RESN,
                "number"_a=first_atom.RESI,
                "index"_a=residues_dict.size()
            ),
            
            "structure"_a=py::dict(
                "atoms"_a=atom_ids,
                "center"_a=py::make_tuple(center[0], center[1], center[2]),
                "n_atoms"_a=atom_ids.size()
            ),
            
            "surface"_a=py::dict(
                "total_sasa"_a=total_sasa,
                "total_area"_a=total_area,
                "standard_sasa"_a=standard_sasa,
                "dsasa"_a=std::max(0.0f, standard_sasa - total_sasa),
                
                "backbone"_a=py::dict(
                    "total"_a=total_components.bb_asa,
                    "polar"_a=total_components.polar_bb_asa,
                    "hydrophobic"_a=total_components.hyd_bb_asa
                ),
                "sidechain"_a=py::dict(
                    "total"_a=total_components.sc_asa,
                    "polar"_a=total_components.polar_sc_asa,
                    "hydrophobic"_a=total_components.hyd_sc_asa
                ),
                
                "groove"_a=py::dict(
                    "major"_a=py::dict(
                        "total"_a=total_components.majorgroove_asa,
                        "polar"_a=total_components.polar_majorgroove_asa,
                        "hydrophobic"_a=total_components.hyd_majorgroove_asa
                    ),
                    "minor"_a=py::dict(
                        "total"_a=total_components.minorgroove_asa,
                        "polar"_a=total_components.polar_minorgroove_asa,
                        "hydrophobic"_a=total_components.hyd_minorgroove_asa
                    ),
                    "none"_a=py::dict(
                        "total"_a=total_components.nogroove_asa,
                        "polar"_a=total_components.polar_nogroove_asa,
                        "hydrophobic"_a=total_components.hyd_nogroove_asa
                    )
                )
            )
        );
        
        // Update chain mapping
        chain_to_residues[first_atom.CHAIN].push_back(residue_id);
    }
    
    // Build chain information
    for (const auto& [chain_id, residue_ids] : chain_to_residues) {
        float total_sasa = 0.0f;
        for (const auto& res_id : residue_ids) {
            auto res = residues_dict[py::str(res_id)];
            total_sasa += res["surface"]["total_sasa"].cast<float>();
        }
        
        // Get chain type from first atom
        std::string chain_type;
        if (!chain_to_atoms[chain_id].empty()) {
            const auto& first_atom_id = chain_to_atoms[chain_id][0];
            const auto& first_atom = atoms[std::stoi(first_atom_id) - 1];
            chain_type = first_atom.MOL_TYPE;
        }
        
        chains_dict[py::str(chain_id)] = py::dict(
            "type"_a=chain_type,
            "residues"_a=residue_ids,
            "surface"_a=py::dict(
                "total_sasa"_a=total_sasa,
                "buried_area"_a=0.0f  // Will be calculated if needed
            )
        );
    }
    
    // Build lookup tables
    for (const auto& [chain_id, atom_ids] : chain_to_atoms) {
        atom_lookup["by_chain"][py::str(chain_id)] = atom_ids;
    }
    
    for (const auto& [residue_id, atom_ids] : residue_to_atoms) {
        atom_lookup["by_residue"][py::str(residue_id)] = atom_ids;
    }
    
    for (size_t i = 0; i < atoms.size(); ++i) {
        if (!atoms[i].ACTIVE) continue;
        atom_lookup["by_id"][py::str(std::to_string(atoms[i].ID))] = i;
    }
    
    for (const auto& [chain_id, residue_ids] : chain_to_residues) {
        res_lookup["by_chain"][py::str(chain_id)] = residue_ids;
    }

    // Create metadata
    std::vector<std::string> chain_ids;
    for (const auto& [chain_id, _] : chain_to_residues) {
        chain_ids.push_back(chain_id);
    }
    
    py::dict metadata = py::dict(
        "total_atoms"_a=atoms.size(),
        "total_residues"_a=residues_dict.size(),
        "chains"_a=chain_ids,
        "types"_a=mol_types,
        "probe_radius"_a=DEFAULT_PROBE_RADIUS
    );
    
    // Assemble final results
    results["atoms"] = atoms_dict;
    results["residues"] = residues_dict;
    results["chains"] = chains_dict;
    results["lookup"] = py::dict(
        "atoms"_a=atom_lookup,
        "residues"_a=res_lookup
    );
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