import dr_sasa_python as sasa
from dr_sasa_python.utils.structure_parser import StructureData, parse_pdb_file, superimpose_structures
import numpy as np
from typing import List, Optional, Dict, Any

def af_structure_to_atom_structs(
    atom_positions: np.ndarray,  # [num_res, num_atom_type, 3]
    atom_mask: np.ndarray,       # [num_res, num_atom_type]
    aatype: np.ndarray,         # [num_res]
    residue_index: np.ndarray,  # [num_res]
    b_factors: Optional[np.ndarray] = None,  # [num_res, num_atom_type]
    chain_id: str = 'A',
    atom_types: List[str] = None,
) -> List[sasa.AtomStruct]:

    # Define standard protein atoms
    atom_types = residue_constants.atom_types

    # Map residue types to three letter codes
    residue_names = dict(enumerate(residue_constants.resnames))

    # Set default b_factors if not provided
    if b_factors is None:
        b_factors = np.zeros_like(atom_mask)

    atom_structs = []
    atom_id = 1

    # Iterate through residues
    for res_idx in range(len(aatype)):
        # Skip residues with no atoms
        if np.sum(atom_mask[res_idx]) < 0.5:
            continue

        resname = residue_names.get(aatype[res_idx], "UNK")

        # Create atom_struct for each atom in residue
        for atom_idx, atom_type in enumerate(atom_types):
            # Skip missing atoms
            if atom_mask[res_idx, atom_idx] < 0.5:
                continue

            coords = atom_positions[res_idx, atom_idx]
            b_factor = b_factors[res_idx, atom_idx]

            # Create atom_struct
            atom_struct = sasa.AtomStruct(
                id=atom_id,
                resi=int(residue_index[res_idx]),
                icode="",  # No insertion codes
                name=atom_type,
                resn=resname,
                chain=chain_id,
                element=atom_type[0],  # First letter of atom type
                structure="protein",
                mol_type="PROTEIN",
                x=float(coords[0]),
                y=float(coords[1]),
                z=float(coords[2]),
                altloc="",
                occupancy=1.0,
                tfactor=float(b_factor),
                charge=""
            )
            atom_structs.append(atom_struct)
            atom_id += 1

    return atom_structs

# Create dummy protein data
num_residues = 5
num_atoms = 37  # Number of standard atoms

# Generate example data
atom_positions = np.random.rand(num_residues, num_atoms, 3)
atom_mask = np.ones((num_residues, num_atoms))
aatype = np.random.randint(0, 20, size=num_residues)
residue_index = np.arange(1, num_residues + 1)
b_factors = np.zeros((num_residues, num_atoms))

# Convert to atom_structs
atom_structs = af_structure_to_atom_structs(
    atom_positions=atom_positions,
    atom_mask=atom_mask,
    aatype=aatype,
    residue_index=residue_index,
    b_factors=b_factors
)
