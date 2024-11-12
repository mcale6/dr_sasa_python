import pytest
import os
import sys
from pathlib import Path
import numpy as np
import dr_sasa_py as sasa
from structure_parser import *
#from bindings.python.utils.structure_parser import StructureData, parse_pdb_file, superimpose_structures

# Test data
@pytest.fixture
def sample_atom_data():
    return {
        'atom_positions': np.array([  # Changed from 'positions' to 'atom_positions'
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [3.0, 0.0, 0.0]
        ]),
        'atom_names': np.array(['N', 'CA', 'C', 'O']),
        'residue_names': np.array(['ALA', 'ALA', 'ALA', 'ALA']),
        'residue_numbers': np.array([1, 1, 1, 1]),
        'chain_ids': np.array(['A', 'A', 'A', 'A']),
        'elements': np.array(['N', 'C', 'C', 'O']),
        'occupancies': np.ones(4),
        'b_factors': np.zeros(4),
        'atom_masks': np.ones(4),
        'alt_locs': np.array([''] * 4),
        'insertion_codes': np.array([''] * 4),
        'charges': np.array([''] * 4)
    }

def create_test_atom():
    """Helper function to create a test atom"""
    return sasa.AtomStruct(
        id=1,
        resi=1,
        icode="",
        name="CA",
        resn="ALA",
        chain="A",
        element="C",
        structure="",
        mol_type="PROTEIN",
        x=0.0,
        y=0.0,
        z=0.0,
        altloc="",
        occupancy=1.0,
        tfactor=0.0,
        charge=""
    )

def test_atom_struct():
    atom = create_test_atom()
    assert atom.ID == 1
    assert atom.NAME == "CA"
    assert atom.RESN == "ALA"
    assert atom.CHAIN == "A"
    assert atom.RESI == 1
    assert list(atom.COORDS) == [0.0, 0.0, 0.0]  # Convert tuple to list for comparison

def test_atom_struct_coordinates():
    atom = create_test_atom()
    # Test setting coordinates
    atom.COORDS = [1.0, 2.0, 3.0]
    assert list(atom.COORDS) == [1.0, 2.0, 3.0]  # Convert tuple to list for comparison
    
    # Test error on wrong number of coordinates
    with pytest.raises(RuntimeError):
        atom.COORDS = [1.0, 2.0]  # Should fail - needs 3 coordinates

def test_structure_data_creation(sample_atom_data):
    structure = StructureData(**sample_atom_data, structure_id="test")
    assert len(structure.atom_names) == 4
    assert structure.structure_id == "test"
    assert np.all(structure.chain_ids == 'A')

def test_structure_validation():
    # Test with invalid coordinates
    bad_data = {
        'atom_positions': np.array([[np.nan, 0.0, 0.0]]),  # Changed from 'positions'
        'atom_names': np.array(['CA']),
        'residue_names': np.array(['ALA']),
        'residue_numbers': np.array([1]),
        'chain_ids': np.array(['A']),
        'elements': np.array(['C']),
        'occupancies': np.ones(1),
        'b_factors': np.zeros(1),
        'atom_masks': np.ones(1),
        'alt_locs': np.array([''] * 1),
        'insertion_codes': np.array([''] * 1),
        'charges': np.array([''] * 1)
    }
    
    with pytest.warns(UserWarning):
        structure = StructureData(**bad_data)
        assert np.all(np.isfinite(structure.atom_positions))

def test_coordinate_transformation(sample_atom_data):
    structure = StructureData(**sample_atom_data)
    
    # Test centering
    centered = structure.transform_coordinates(center=True)
    assert np.allclose(np.mean(centered.atom_positions, axis=0), 0.0)
    
    # Test rotation
    rotation = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]])  # 90° around z
    rotated = structure.transform_coordinates(rotation=rotation)
    assert np.allclose(rotated.atom_positions[1], [0.0, 1.0, 0.0])

def test_atom_selection(sample_atom_data):
    structure = StructureData(**sample_atom_data)
    
    # Select CA atoms
    ca_mask = structure.atom_names == 'CA'
    ca_struct = structure.select_atoms(ca_mask)
    assert len(ca_struct.atom_names) == 1
    assert ca_struct.atom_names[0] == 'CA'

# Existing tests
def test_pdb_parsing():
    if not os.path.exists("tests/data/3i40.pdb"):
        pytest.skip("Test PDB file not found")
    structure = parse_pdb_file("tests/data/3i40.pdb")
    assert len(structure.atom_names) > 0
    assert np.all(np.isfinite(structure.atom_positions))
    assert np.all(structure.occupancies > 0)

def test_simple_sasa_creation():
    calc = sasa.SimpleSASA()
    assert calc is not None
    
    calc_custom = sasa.SimpleSASA(probe_radius=1.2)
    assert calc_custom is not None

@pytest.mark.skipif(not os.path.exists("tests/data/3i40.pdb"), reason="Test PDB file not found")
def test_simple_sasa_calculation():
    calc = sasa.SimpleSASA()
    results = calc.calculate("tests/data/3i40.pdb")
    assert isinstance(results, dict)
    assert len(results) > 0

def test_structure_superposition(sample_atom_data):
    structure = StructureData(**sample_atom_data)
    
    # Add small random displacement to create second structure
    displaced = structure.transform_coordinates(
        translation=np.array([1.0, 2.0, -1.0]),
        rotation=np.eye(3)
    )
    
    # Superimpose
    aligned = superimpose_structures(displaced, structure)
    
    # Check RMSD is small
    rmsd = np.sqrt(np.mean(np.sum((aligned.atom_positions - 
                                  structure.atom_positions)**2, axis=1)))
    assert rmsd < 1e-10

if __name__ == "__main__":
    pytest.main([__file__, "-v"])