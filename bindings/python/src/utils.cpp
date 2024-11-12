#include "utils.hpp"

static const std::map<std::string, float> STANDARD_SASA_VALUES = {
    {"ALA", 92.40211452},
    {"GLN", 184.71688},
    {"TRP", 235.3483229},
    {"SER", 122.353095},
    {"ARG", 219.941475},
    {"ILE", 158.07277},
    {"ASN", 146.06073},
    {"ASP", 154.71124},
    {"HIS", 194.617376},
    {"MET", 191.547244},
    {"LYS", 201.689792},
    {"LEU", 169.04496452},
    {"THR", 145.526463},
    {"PHE", 201.7065},
    {"TYR", 209.07715},
    {"GLU", 160.482561},
    {"CYS", 130.28853},
    {"PRO", 126.877028},
    {"GLY", 80.23533},
    {"VAL", 138.90233},
    {"A", 163.5912717},
    {"C", 162.86302775},
    {"T", 173.03036075},
    {"G", 164.253836},
    {"DA", 163.5912717},
    {"DC", 162.86302775},
    {"DT", 173.03036075},
    {"DG", 164.253836}
};

void RelativeSASA2(std::vector<atom_struct>& atoms) {
    for (auto& atom : atoms) {
        if (!atom.ACTIVE) continue;
        
        auto it = STANDARD_SASA_VALUES.find(atom.RESN);
        if (it != STANDARD_SASA_VALUES.end()) {
            atom.SASA /= it->second;
        } else {
            atom.SASA = 0.0f;  // ? or another default value for unknown residues
        }
    }
}


py::dict create_analysis_results_(const std::vector<atom_struct>& atoms, bool include_matrix) {
    py::dict results;
    using shape_t = std::vector<ssize_t>;
    
    // Process atom-level data (this part remains serial as it's creating Python objects)
    for (const auto& atom : atoms) {
        if (!atom.ACTIVE) continue;

        py::dict atom_data = py::dict(
            "name"_a=atom.NAME,
            "resname"_a=atom.RESN,
            "chain"_a=atom.CHAIN,
            "resid"_a=atom.RESI,
            "struct_type"_a=atom.STRUCT_TYPE,
            "coords"_a=py::make_tuple(atom.COORDS[0], atom.COORDS[1], atom.COORDS[2]),
            "sphere_area"_a=atom.AREA,
            "sasa"_a=atom.SASA,
            "polar"_a=atom.POLAR,
            "charge"_a=atom.CHARGE
        );

        // Contacts and overlaps processing...
        // [existing code]

        results[py::str(std::to_string(atom.ID))] = std::move(atom_data);
    }

    // Prepare residue data
    std::map<std::tuple<std::string, std::string, int>, size_t> residue_indices;
    size_t n_residues = 0;
    
    // First pass to count residues (needs to be serial)
    for (const auto& atom : atoms) {
        if (!atom.ACTIVE) continue;
        auto res_key = std::make_tuple(atom.CHAIN, atom.RESN, atom.RESI);
        if (residue_indices.find(res_key) == residue_indices.end()) {
            residue_indices[res_key] = n_residues++;
        }
    }

    // Calculate relative SASA
    auto atoms_copy = atoms;
    RelativeSASA2(atoms_copy);

    // Allocate result arrays
    std::vector<float> sasa_values(n_residues, 0.0f);
    std::vector<float> area_values(n_residues, 0.0f);
    std::vector<float> rel_sasa_values(n_residues, 0.0f);
    std::vector<float> coordinates(n_residues * 3, 0.0f);
    std::vector<std::string> chains(n_residues);
    std::vector<std::string> resnames(n_residues);
    std::vector<int> resids(n_residues);
    std::vector<int> atom_counts(n_residues, 0);

    // Parallel processing of residue data
    #pragma omp parallel
    {
        // Thread-local storage
        //const int thread_id = omp_get_thread_num();
        //const int n_threads = omp_get_num_threads();
        
        // Each thread gets its own accumulation arrays
        std::vector<float> local_sasa(n_residues, 0.0f);
        std::vector<float> local_area(n_residues, 0.0f);
        std::vector<float> local_rel_sasa(n_residues, 0.0f);
        std::vector<float> local_coords(n_residues * 3, 0.0f);
        std::vector<int> local_counts(n_residues, 0);

        // Process atoms in parallel
        #pragma omp for schedule(guided)
        for (size_t i = 0; i < atoms.size(); ++i) {
            const auto& atom = atoms[i];
            if (!atom.ACTIVE) continue;

            auto res_key = std::make_tuple(atom.CHAIN, atom.RESN, atom.RESI);
            size_t res_idx = residue_indices[res_key];
            
            // Accumulate in thread-local storage
            local_sasa[res_idx] += atom.SASA;
            local_area[res_idx] += atom.AREA;
            local_rel_sasa[res_idx] += atoms_copy[i].SASA;
            
            local_coords[res_idx * 3] += atom.COORDS[0];
            local_coords[res_idx * 3 + 1] += atom.COORDS[1];
            local_coords[res_idx * 3 + 2] += atom.COORDS[2];
            
            local_counts[res_idx]++;

            // First thread to process a residue sets its metadata
            if (local_counts[res_idx] == 1) {
                #pragma omp critical(metadata)
                {
                    if (atom_counts[res_idx] == 0) {
                        chains[res_idx] = atom.CHAIN;
                        resnames[res_idx] = atom.RESN;
                        resids[res_idx] = atom.RESI;
                    }
                }
            }
        }

        // Combine results from all threads
        #pragma omp critical(accumulate)
        {
            for (size_t i = 0; i < n_residues; ++i) {
                sasa_values[i] += local_sasa[i];
                area_values[i] += local_area[i];
                rel_sasa_values[i] += local_rel_sasa[i];
                coordinates[i * 3] += local_coords[i * 3];
                coordinates[i * 3 + 1] += local_coords[i * 3 + 1];
                coordinates[i * 3 + 2] += local_coords[i * 3 + 2];
                atom_counts[i] += local_counts[i];
            }
        }
    }

    // Calculate averages in parallel
    #pragma omp parallel for schedule(static)
    for (size_t i = 0; i < n_residues; ++i) {
        if (atom_counts[i] > 0) {
            float inv_count = 1.0f / atom_counts[i];
            rel_sasa_values[i] *= inv_count;
            coordinates[i * 3] *= inv_count;
            coordinates[i * 3 + 1] *= inv_count;
            coordinates[i * 3 + 2] *= inv_count;
        }
    }

    // Create residue data with numpy arrays (serial)
    py::dict residue_data;
    residue_data["sasa"] = py::array_t<float>(shape_t{(ssize_t)n_residues}, sasa_values.data());
    residue_data["area"] = py::array_t<float>(shape_t{(ssize_t)n_residues}, area_values.data());
    residue_data["relative_sasa"] = py::array_t<float>(shape_t{(ssize_t)n_residues}, rel_sasa_values.data());
    residue_data["coordinates"] = py::array_t<float>(shape_t{(ssize_t)n_residues, 3}, coordinates.data());
    residue_data["chains"] = py::cast(chains);
    residue_data["resnames"] = py::cast(resnames);
    residue_data["resids"] = py::array_t<int>(shape_t{(ssize_t)n_residues}, resids.data());
    residue_data["n_atoms"] = py::array_t<int>(shape_t{(ssize_t)n_residues}, atom_counts.data());
    
    results["residue_data"] = residue_data;
    
    return results;
}


py::dict create_analysis_results(const std::vector<atom_struct>& atoms, bool include_matrix) {
    py::dict results;
    
    // First collect per-residue data
    std::map<std::tuple<std::string, std::string, int>, py::dict> residue_data;
    
    // Calculate relative SASA for all atoms
    auto atoms_copy = atoms;
    RelativeSASA2(atoms_copy);
    
    // Process atoms and aggregate residue data
    for (size_t i = 0; i < atoms.size(); ++i) {
        const auto& atom = atoms[i];
        if (!atom.ACTIVE) continue;

        // Create atom data as before
        py::dict atom_data = py::dict(
            "name"_a=atom.NAME,
            "resname"_a=atom.RESN,
            "chain"_a=atom.CHAIN,
            "resid"_a=atom.RESI,
            "struct_type"_a=atom.STRUCT_TYPE,
            "coords"_a=py::make_tuple(atom.COORDS[0], atom.COORDS[1], atom.COORDS[2]),
            "sphere_area"_a=atom.AREA,
            "sasa"_a=atom.SASA,
            "polar"_a=atom.POLAR,
            "charge"_a=atom.CHARGE,
            "relative_sasa"_a=atoms_copy[i].SASA  // Get relative SASA from copied atoms
        );

        // Add contacts data
        py::dict contacts_data;
        for (uint32_t pos : atom.INTERACTION_SASA_P) {
            const auto& other_atom = atoms[pos];
            contacts_data[py::str(std::to_string(other_atom.ID))] = py::dict(
                "struct_type"_a=other_atom.STRUCT_TYPE,
                "contact_area"_a=atom.CONTACT_AREA.at(other_atom.ID),
                "distance"_a=atom.DISTANCES.at(pos)
            );
        }
        atom_data["contacts"] = contacts_data;

        // Add overlaps data
        std::vector<py::dict> overlaps_data;
        for (size_t i = 0; i < atom.ov_table.size(); ++i) {
            const auto& overlap_set = atom.ov_table[i];
            std::vector<uint32_t> overlap_ids;
            for (uint32_t pos : overlap_set) {
                overlap_ids.push_back(atoms[pos].ID);
            }
            
            overlaps_data.push_back(py::dict(
                "atoms"_a=overlap_ids,
                "overlap_area"_a=atom.ov_table_area[i],
                "normalized_area"_a=atom.ov_norm_area[i],
                "buried_area"_a=atom.AREA_BURIED_BY_ATOM_area[i]
            ));
        }
        atom_data["overlaps"] = overlaps_data;

        results[py::str(std::to_string(atom.ID))] = std::move(atom_data);

        // Aggregate residue data
        auto res_key = std::make_tuple(atom.CHAIN, atom.RESN, atom.RESI);
        if (residue_data.find(res_key) == residue_data.end()) {
            residue_data[res_key] = py::dict(
                "chain"_a=atom.CHAIN,
                "resname"_a=atom.RESN,
                "resid"_a=atom.RESI,
                "total_sasa"_a=atom.SASA,
                "total_area"_a=atom.AREA,
                "relative_sasa"_a=atoms_copy[i].SASA,
                "n_atoms"_a=1,
                "center"_a=py::make_tuple(atom.COORDS[0], atom.COORDS[1], atom.COORDS[2])
            );
        } else {
            auto& res = residue_data[res_key];
            res["total_sasa"] = res["total_sasa"].cast<float>() + atom.SASA;
            res["total_area"] = res["total_area"].cast<float>() + atom.AREA;
            res["relative_sasa"] = res["relative_sasa"].cast<float>() + atoms_copy[i].SASA;
            res["n_atoms"] = res["n_atoms"].cast<int>() + 1;
            
            // Update center
            auto center = res["center"].cast<py::tuple>();
            res["center"] = py::make_tuple(
                center[0].cast<float>() + atom.COORDS[0],
                center[1].cast<float>() + atom.COORDS[1],
                center[2].cast<float>() + atom.COORDS[2]
            );
        }
    }

    // Finalize residue data (calculate averages)
    py::list res_data_list;
    for (auto& [key, res] : residue_data) {
        int n_atoms = res["n_atoms"].cast<int>();
        auto center = res["center"].cast<py::tuple>();
        
        // Calculate average center
        res["center"] = py::make_tuple(
            center[0].cast<float>() / n_atoms,
            center[1].cast<float>() / n_atoms,
            center[2].cast<float>() / n_atoms
        );
        // Calculate average relative SASA
        res["relative_sasa"] = res["relative_sasa"].cast<float>() / n_atoms;
        res_data_list.append(res);
    }

    results["residue_data"] = res_data_list;
    
    return results;
}

py::dict generate_interaction_matrices(std::vector<atom_struct>& atoms) {
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

py::dict generate_intra_matrices(std::vector<atom_struct>& atoms) {
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