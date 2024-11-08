// bindings/python/src/utils.cpp
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