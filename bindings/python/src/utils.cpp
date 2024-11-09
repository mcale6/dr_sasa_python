#include "utils.hpp"

py::dict create_analysis_results(const std::vector<atom_struct>& atoms, bool include_matrix) {
    py::dict results;
    
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
    }
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