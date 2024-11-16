#include "common.hpp"

static const std::map<std::string, float> STANDARD_SASA_VALUES = {
                            {"ALA",92.40211452},
                            {"GLN",184.71688},
                            {"TRP",235.3483229},
                            {"SER",122.353095},
                            {"ARG",219.941475},
                            {"ILE",158.07277},
                            {"ASN",146.06073},
                            {"ASP",154.71124},
                            {"HIS",194.617376},
                            {"MET",191.547244},
                            {"LYS",201.689792},
                            {"LEU",169.04496452},
                            {"THR",145.526463}, 
                            {"PHE",201.7065},
                            {"TYR",209.07715},
                            {"GLU",160.482561},
                            {"CYS",130.28853},
                            {"PRO",126.877028},
                            {"GLY",80.23533}, 
                            {"VAL",138.90233},  
                            { "A",163.5912717},
                            { "C",162.86302775},
                            { "T",173.03036075},
                            { "G",164.253836},
                            { "DA",163.5912717},
                            { "DC",162.86302775},
                            { "DT",173.03036075},
                            { "DG",164.253836}
};
py::dict create_analysis_results(const std::vector<atom_struct>& atoms, bool include_matrix) {
    py::dict results;
    py::dict atom_data;
    py::list residue_data;  // Keep as list
    py::dict residue_index; // New index structure
    
    // collect per-residue data with contacts and overlaps
    std::map<std::tuple<std::string, std::string, int>, py::dict> residue_temp_data;
    
    // Process atoms and aggregate residue data
    for (size_t i = 0; i < atoms.size(); ++i) {
        const auto& atom = atoms[i];
        if (!atom.ACTIVE) continue;

        // Create atom data with new structure
        py::dict atom_info = py::dict(
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