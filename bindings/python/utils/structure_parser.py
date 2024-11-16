import dataclasses
from pathlib import Path
import warnings
from typing import Dict, List, Optional, Tuple, Union
import numpy as np
import sys
import dr_sasa_py as sasa

# PDB standard chain IDs
PDB_CHAIN_IDS = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789'

# Standard residue mapping
RESIDUE_NAMES_3TO1 = {
    'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
    'GLN': 'Q', 'GLU': 'E', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
    'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'PHE': 'F', 'PRO': 'P',
    'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
    'UNK': 'X'
}

@dataclasses.dataclass
class StructureData:
    """Structure representation with validation and transformation capabilities."""
    # Required fields
    atom_positions: np.ndarray  # [num_atoms, 3] - Cartesian coordinates of atoms
    atom_names: np.ndarray     # [num_atoms] - PDB atom names (e.g., 'CA', 'N')
    residue_names: np.ndarray  # [num_atoms] - Three letter residue codes
    chain_ids: np.ndarray      # [num_atoms] - Chain identifiers
    residue_numbers: np.ndarray # [num_atoms] - Residue sequence numbers
    elements: np.ndarray       # [num_atoms] - Element symbols
    
    # Optional fields with defaults
    occupancies: np.ndarray = dataclasses.field(default_factory=lambda: np.ones(0))      # [num_atoms] - Occupancy values
    b_factors: np.ndarray = dataclasses.field(default_factory=lambda: np.zeros(0))       # [num_atoms] - Temperature factors
    atom_masks: np.ndarray = dataclasses.field(default_factory=lambda: np.ones(0))       # [num_atoms] - Mask for each atom
    alt_locs: np.ndarray = dataclasses.field(default_factory=lambda: np.array([]))       # [num_atoms] - Alternate location indicators
    insertion_codes: np.ndarray = dataclasses.field(default_factory=lambda: np.array([])) # [num_atoms] - PDB insertion codes
    charges: np.ndarray = dataclasses.field(default_factory=lambda: np.array([]))        # [num_atoms] - Atomic charges
    structure_id: str = "structure"  # Identifier for the structure

    def __post_init__(self):
        """Validate and initialize structure data after creation."""
        n_atoms = len(self.atom_names)
        
        # Initialize any empty arrays with correct size
        if len(self.occupancies) == 0:
            self.occupancies = np.ones(n_atoms)
        if len(self.b_factors) == 0:
            self.b_factors = np.zeros(n_atoms)
        if len(self.atom_masks) == 0:
            self.atom_masks = np.ones(n_atoms)
        if len(self.alt_locs) == 0:
            self.alt_locs = np.array([''] * n_atoms)
        if len(self.insertion_codes) == 0:
            self.insertion_codes = np.array([''] * n_atoms)
        if len(self.charges) == 0:
            self.charges = np.array([''] * n_atoms)
            
        self.validate()

    def validate(self) -> None:
        """Validate structure data and fix common issues."""
        n_atoms = len(self.atom_names)
        
        # Check array shapes
        if self.atom_positions.shape != (n_atoms, 3):
            raise ValueError(f"atom_positions should have shape ({n_atoms}, 3), got {self.atom_positions.shape}")
        
        # Check array lengths
        for name, arr in [
            ('atom_names', self.atom_names),
            ('residue_names', self.residue_names),
            ('chain_ids', self.chain_ids),
            ('residue_numbers', self.residue_numbers),
            ('elements', self.elements),
            ('occupancies', self.occupancies),
            ('b_factors', self.b_factors),
            ('atom_masks', self.atom_masks),
            ('alt_locs', self.alt_locs),
            ('insertion_codes', self.insertion_codes),
            ('charges', self.charges)
        ]:
            if len(arr) != n_atoms:
                raise ValueError(f"{name} should have length {n_atoms}, got {len(arr)}")

        # Fix non-finite coordinates
        if np.any(~np.isfinite(self.atom_positions)):
            warnings.warn("Found non-finite coordinates, setting to 0.0")
            self.atom_positions[~np.isfinite(self.atom_positions)] = 0.0

        # Validate and fix chain IDs
        invalid_chains = ~np.isin(self.chain_ids, list(PDB_CHAIN_IDS))
        if np.any(invalid_chains):
            warnings.warn(f"Found invalid chain IDs, replacing with 'X': {set(self.chain_ids[invalid_chains])}")
            self.chain_ids[invalid_chains] = 'X'

    @classmethod
    def from_arrays(cls, atom_data: Dict[str, np.ndarray], structure_id: str = "array_structure") -> 'StructureData':
        """Create StructureData from dictionary of numpy arrays."""
        # Handle coordinates
        positions = atom_data.get('atom_positions')
        if positions is None:
            # Try to construct from x, y, z arrays
            x = atom_data.get('x', np.zeros(len(atom_data['atom_names'])))
            y = atom_data.get('y', np.zeros(len(atom_data['atom_names'])))
            z = atom_data.get('z', np.zeros(len(atom_data['atom_names'])))
            positions = np.stack([x, y, z], axis=-1)
            
        n_atoms = len(atom_data['atom_names'])
        
        # Prepare kwargs with required and optional fields
        kwargs = {
            'atom_positions': positions,
            'atom_names': atom_data['atom_names'],
            'residue_names': atom_data['residue_names'],
            'chain_ids': atom_data.get('chain_ids', np.array(['A'] * n_atoms)),
            'residue_numbers': atom_data['residue_numbers'],
            'elements': atom_data.get('elements', np.array([name[0] for name in atom_data['atom_names']])),
            'structure_id': structure_id
        }
        
        # Add optional fields if present
        for field in ['occupancies', 'b_factors', 'atom_masks', 'alt_locs', 'insertion_codes', 'charges']:
            if field in atom_data:
                kwargs[field] = atom_data[field]
                
        return cls(**kwargs)

    def transform_coordinates(self, 
                            rotation: Optional[np.ndarray] = None,
                            translation: Optional[np.ndarray] = None,
                            center: bool = False) -> 'StructureData':
        """Apply coordinate transformations."""
        new_positions = self.atom_positions.copy()
        
        if center:
            centroid = np.mean(new_positions, axis=0)
            new_positions -= centroid
            
        if rotation is not None:
            if rotation.shape != (3, 3):
                raise ValueError("Rotation matrix must be 3x3")
            new_positions = new_positions @ rotation.T
            
        if translation is not None:
            translation = np.asarray(translation).reshape(1, 3)
            new_positions += translation
            
        new_data = dataclasses.replace(self, atom_positions=new_positions)
        return new_data

    def select_atoms(self, mask: np.ndarray) -> 'StructureData':
        """Create new StructureData with only selected atoms."""
        if len(mask) != len(self.atom_names):
            raise ValueError("Selection mask must match number of atoms")
        
        # Create new structure with selected atoms
        return dataclasses.replace(
            self,
            atom_positions=self.atom_positions[mask],
            atom_names=self.atom_names[mask],
            residue_names=self.residue_names[mask],
            chain_ids=self.chain_ids[mask],
            residue_numbers=self.residue_numbers[mask],
            elements=self.elements[mask],
            occupancies=self.occupancies[mask],
            b_factors=self.b_factors[mask],
            atom_masks=self.atom_masks[mask],
            alt_locs=self.alt_locs[mask],
            insertion_codes=self.insertion_codes[mask],
            charges=self.charges[mask]
        )

    def to_atom_structs(self) -> List[sasa.AtomStruct]:
        """Convert to list of AtomStruct objects."""
        atoms = []
        for i in range(len(self.atom_names)):
            if self.atom_masks[i] < 0.5:
                continue

            atom = sasa.AtomStruct(
                id=i + 1,
                resi=int(self.residue_numbers[i]),
                icode=str(self.insertion_codes[i]),
                name=str(self.atom_names[i]),
                resn=str(self.residue_names[i]),
                chain=str(self.chain_ids[i]),
                element=str(self.elements[i]),
                structure=self.structure_id,
                mol_type="",  # Will be determined later
                x=float(self.atom_positions[i, 0]),
                y=float(self.atom_positions[i, 1]),
                z=float(self.atom_positions[i, 2]),
                altloc=str(self.alt_locs[i]),
                occupancy=float(self.occupancies[i]),
                tfactor=float(self.b_factors[i]),
                charge=str(self.charges[i])
            )
            atoms.append(atom)
        return atoms

def superimpose_structures(mobile: StructureData, 
                         target: StructureData,
                         atom_mask: Optional[np.ndarray] = None) -> StructureData:
    """Superimpose mobile structure onto target structure."""
    if atom_mask is None:
        atom_mask = np.ones(len(mobile.atom_names), dtype=bool)
        
    if len(mobile.atom_positions) != len(target.atom_positions):
        raise ValueError("Structures must have same number of atoms")
        
    # Get coordinates for superposition
    mobile_coords = mobile.atom_positions[atom_mask]
    target_coords = target.atom_positions[atom_mask]
    
    # Center the structures
    mobile_center = np.mean(mobile_coords, axis=0)
    target_center = np.mean(target_coords, axis=0)
    
    mobile_centered = mobile_coords - mobile_center
    target_centered = target_coords - target_center
    
    # Calculate rotation matrix using SVD
    covariance = mobile_centered.T @ target_centered
    U, _, Vt = np.linalg.svd(covariance)
    rotation = Vt.T @ U.T
    
    # Fix chirality if needed
    if np.linalg.det(rotation) < 0:
        Vt[-1] *= -1
        rotation = Vt.T @ U.T
    
    # Apply transformation
    translation = target_center - mobile_center @ rotation
    
    return mobile.transform_coordinates(rotation=rotation, translation=translation)


def parse_pdb_file(file_path: Union[str, Path]) -> StructureData:
    """Parse PDB file into StructureData."""
    file_path = Path(file_path)
    if not file_path.exists():
        raise FileNotFoundError(f"PDB file not found: {file_path}")
        
    # Initialize lists to store atom data
    data = {
        'atom_names': [], 'residue_names': [], 'chain_ids': [],
        'residue_numbers': [], 'elements': [], 'occupancies': [],
        'b_factors': [], 'alt_locs': [], 'insertion_codes': [],
        'charges': [], 'positions': []
    }
    
    with open(file_path) as f:
        for line in f:
            if line.startswith(('ATOM', 'HETATM')):
                try:
                    # Parse PDB line following standard format
                    data['atom_names'].append(line[12:16].strip())
                    data['residue_names'].append(line[17:20].strip())
                    data['chain_ids'].append(line[21:22].strip() or 'A')
                    data['residue_numbers'].append(int(line[22:26]))
                    data['alt_locs'].append(line[16:17].strip())
                    data['insertion_codes'].append(line[26:27].strip())
                    
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    data['positions'].append([x, y, z])
                    
                    data['occupancies'].append(float(line[54:60] or '1.0'))
                    data['b_factors'].append(float(line[60:66] or '0.0'))
                    data['elements'].append(line[76:78].strip() or line[12:16].strip()[0])
                    data['charges'].append(line[78:80].strip())
                    
                except Exception as e:
                    warnings.warn(f"Failed to parse line: {line.strip()}\nError: {e}")
                    continue

    # Convert lists to numpy arrays
    arrays = {
        'atom_positions': np.array(data['positions']),
        'atom_names': np.array(data['atom_names']),
        'residue_names': np.array(data['residue_names']),
        'chain_ids': np.array(data['chain_ids']),
        'residue_numbers': np.array(data['residue_numbers']),
        'elements': np.array(data['elements']),
        'occupancies': np.array(data['occupancies']),
        'b_factors': np.array(data['b_factors']),
        'atom_masks': np.ones(len(data['atom_names'])),
        'alt_locs': np.array(data['alt_locs']),
        'insertion_codes': np.array(data['insertion_codes']),
        'charges': np.array(data['charges']),
        'structure_id': file_path.stem
    }
    
    return StructureData(**arrays)


def create_atom_structs(
    atom_data: Dict[str, np.ndarray],
    structure_id: Optional[str] = None
) -> List[sasa.AtomStruct]:
    """Create AtomStruct list from array data.
    
    Args:
        atom_data: Dictionary containing arrays describing the structure:
            Required:
            - atom_names: [N] atom names
            - residue_names: [N] residue names
            - residue_numbers: [N] residue sequence numbers
            And either:
            - positions: [N, 3] coordinates or
            - x, y, z: [N] coordinate arrays
            
            Optional:
            - chain_ids: [N] chain identifiers (default: 'A')
            - elements: [N] element symbols (default: first letter of atom name) ## improve
            - occupancies: [N] occupancy values (default: 1.0)
            - b_factors: [N] B-factors (default: 0.0)
            - atom_masks: [N] mask for each atom (default: 1)
            - alt_locs: [N] alternate location indicators (default: '')
            - insertion_codes: [N] insertion codes (default: '')
            - charges: [N] atomic charges (default: '')
            
        structure_id: Optional identifier for the structure
        
    Returns:
        List of AtomStruct objects ready for SASA calculation
    """
    struct_data = StructureData.from_arrays(atom_data, structure_id or "structure")
    return struct_data.to_atom_structs()
