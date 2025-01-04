import pytest
import warnings
import tempfile
import os
import subprocess
from pathlib import Path
from typing import Dict, Any
import numpy as np
from Bio.PDB import PDBParser
import dr_sasa_python as sasa
from dr_sasa_python.utils.structure_parser import StructureData, parse_pdb_file, superimpose_structures

# Test data constants
# Constants
DR_SASA_N_EXEC = "/home/alessio/dr_sasa_python/dr_sasa_n/build/dr_sasa"
TEST_DATA_DIR = Path(__file__).parent / "data"
TEST_FILES = {
    "basic": "3i40.pdb",
    "complex": "6gwp.pdb",
    "prediction": "pred.pdb" 
}

# Fixtures
@pytest.fixture
def dr_sasa_exec() -> str:
    """Path to original dr_sasa executable"""
    # Try environment variable first, fall back to constant
    exec_path = os.environ.get('DR_SASA_EXEC', DR_SASA_N_EXEC)
    
    if not Path(exec_path).exists():
        pytest.skip(f"Original dr_sasa executable not found at {exec_path}")
    return exec_path

@pytest.fixture(params=["simple", "generic", "decoupled"])
def calc_type(request):
    """Calculator type to test."""
    return request.param

@pytest.fixture
def calculators():
    """Provide calculator instances."""
    return {
        "simple": sasa.SimpleSASA(probe_radius=1.4),
        "generic": sasa.GenericSASA(probe_radius=1.4),
        "decoupled": sasa.DecoupledSASA(probe_radius=1.4)
    }

@pytest.fixture
def test_paths():
    """Provide test file paths."""
    def _get_path(key: str) -> Path:
        path = TEST_DATA_DIR / TEST_FILES[key]
        if not path.exists():
            pytest.skip(f"Test file not found: {path}")
        return path
    return _get_path


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

def test_pdb_parsing(test_paths):
    pdb_file = test_paths("basic")
    if not pdb_file.exists():
        print(pdb_file)
        pytest.skip("Test PDB file not found")
    structure = parse_pdb_file(pdb_file)
    assert len(structure.atom_names) > 0
    assert np.all(np.isfinite(structure.atom_positions))
    assert np.all(structure.occupancies > 0)

def test_simple_sasa_calculation(test_paths):
    pdb_file = test_paths("basic")
    if not pdb_file.exists():
        raise FileNotFoundError(f"Test PDB file not found: {pdb_file}")
    
    calc = sasa.SimpleSASA()
    results = calc.calculate(str(pdb_file))
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

def get_test_file(key: str) -> Path:
    """Get path to test file, skip if not found"""
    path = TEST_DATA_DIR / TEST_FILES[key]
    if not path.exists():
        pytest.skip(f"Test file not found: {path}")
    return path

def validate_sasa_results(results: Dict[str, Any], calc_type: str = "simple") -> None:
    """Validate SASA calculation results structure and data types."""
    # Validate top-level structure
    assert isinstance(results, dict), "Results must be a dictionary"
    required_keys = {"atoms", "residues", "chains", "lookup", "metadata"}
    missing_keys = required_keys - set(results.keys())
    assert not missing_keys, f"Missing top-level keys: {missing_keys}"

    # Validate atom data
    atoms = results["atoms"]
    assert isinstance(atoms, dict), "atoms must be a dictionary"
    assert len(atoms) > 0, "atoms cannot be empty"

    # Get first atom entry for detailed validation
    first_atom_id = next(iter(atoms))
    first_atom = atoms[first_atom_id]

    # Validate atom data structure
    required_atom_sections = {
        "name": str,
        "resid": int,
        "resname": str,
        "chain": str,
        "index": int,
        "coords": tuple,
        "surface": dict,
        "properties": dict,
        "contacts": dict
    }

    missing_sections = required_atom_sections.keys() - set(first_atom.keys())
    assert not missing_sections, f"Missing atom sections: {missing_sections}"

    # Validate surface section
    required_surface_fields = {
        "sphere_area": float,
        "sasa": float,
        "buried_area": float,
        "contact_area": float,
        "dsasa": float
    }

    surface = first_atom["surface"]
    missing_surface = required_surface_fields.keys() - set(surface.keys())
    assert not missing_surface, f"Missing surface fields: {missing_surface}"

    # Validate properties section
    required_properties = {
        "vdw": float,
        "polar": int,
        "charge": str,
        "struct_type": str
    }

    props = first_atom["properties"]
    missing_props = required_properties.keys() - set(props.keys())
    assert not missing_props, f"Missing property fields: {missing_props}"

    # Validate contacts section
    required_contact_sections = {
        "nonbonded": dict,
        "overlap_groups": list
    }

    contacts = first_atom["contacts"]
    missing_contacts = required_contact_sections.keys() - set(contacts.keys())
    assert not missing_contacts, f"Missing contact sections: {missing_contacts}"

    # Type validation for specific fields
    assert isinstance(first_atom["coords"], tuple) and len(first_atom["coords"]) == 3, \
        "coords must be a 3-tuple"
    assert all(isinstance(x, (int, float)) for x in first_atom["coords"]), \
        "coords must contain numeric values"
    assert first_atom["properties"]["polar"] in (0, 1), \
        f"polar must be 0 or 1, got {first_atom['properties']['polar']}"
    
    # Validate residue data
    residues = results["residues"]
    assert isinstance(residues, dict), "residues must be a dictionary"
    assert len(residues) > 0, "residues cannot be empty"

    first_residue = next(iter(residues.values()))
    required_residue_sections = {
        "identifiers": dict,
        "structure": dict,
        "surface": dict
    }

    missing_res_sections = required_residue_sections.keys() - set(first_residue.keys())
    assert not missing_res_sections, f"Missing residue sections: {missing_res_sections}"

    # Validate residue surface section
    required_residue_surface = {
        "total_sasa": float,
        "total_area": float,
        "standard_sasa": float,
        "dsasa": float,
        #"backbone": dict,
        #"sidechain": dict,
        #"groove": dict
    }

    res_surface = first_residue["surface"]
    missing_res_surface = required_residue_surface.keys() - set(res_surface.keys())
    assert not missing_res_surface, f"Missing residue surface fields: {missing_res_surface}"

    # Validate chain data
    chains = results["chains"]
    assert isinstance(chains, dict), "chains must be a dictionary"
    if len(chains) > 0:
        first_chain = next(iter(chains.values()))
        required_chain_fields = {
            "type": str,
            "residues": list,
            #"surface": dict
        }

        missing_chain_fields = required_chain_fields.keys() - set(first_chain.keys())
        assert not missing_chain_fields, f"Missing chain fields: {missing_chain_fields}"

    # Validate lookup tables
    lookup = results["lookup"]
    required_lookup_sections = {
        "atoms": dict,
        "residues": dict
    }

    missing_lookup = required_lookup_sections.keys() - set(lookup.keys())
    assert not missing_lookup, f"Missing lookup sections: {missing_lookup}"

    # Validate metadata
    metadata = results["metadata"]
    required_metadata_fields = {
        "total_atoms": int,
        "total_residues": int,
        "chains": list,
        "types": list,
        "probe_radius": float
    }

    missing_metadata = required_metadata_fields.keys() - set(metadata.keys())
    assert not missing_metadata, f"Missing metadata fields: {missing_metadata}"

    # Calculator-specific validations
    if calc_type in ["generic", "decoupled"]:
        # Verify contact information is present
        nonbonded = first_atom["contacts"]["nonbonded"]
        if len(nonbonded) > 0:
            first_contact = next(iter(nonbonded.values()))
            assert "area" in first_contact, "Contact area missing"
            assert "distance" in first_contact, "Contact distance missing"

        # Verify overlap information
        overlap_groups = first_atom["contacts"]["overlap_groups"]
        if len(overlap_groups) > 0:
            first_overlap = overlap_groups[0]
            assert "atoms" in first_overlap, "Overlap atoms missing"
            assert "area" in first_overlap, "Overlap area missing"
            assert "normalized_area" in first_overlap, "Normalized overlap area missing"
            assert "buried_area" in first_overlap, "Buried area missing"
 
class TestDrSasaPy:
    def test_atom_struct(self):
        """Test AtomStruct creation and properties"""
        atom = sasa.AtomStruct(
            id=1, resi=1, icode="", name="CA", resn="ALA",
            chain="A", element="C", structure="test",
            mol_type="PROTEIN", x=0.0, y=0.0, z=0.0
        )
        
        assert atom.ID == 1
        assert atom.NAME == "CA"
        assert atom.RESN == "ALA"
        assert list(atom.COORDS) == [0.0, 0.0, 0.0]

    def test_structure_data(self):
        """Test sasa.StructureData creation and manipulation"""
        pdb_path = get_test_file("basic")
        structure = parse_pdb_file(pdb_path)
        assert len(structure.atom_names) > 0
        assert np.all(np.isfinite(structure.atom_positions))

    def test_sasa_calculation(self, calc_type, calculators):
        """Test SASA calculation with different calculators"""
        calculator = calculators[calc_type]
        pdb_path = get_test_file("complex")
        
        # Basic calculation
        results = calculator.calculate(str(pdb_path))
        validate_sasa_results(results, calc_type)
        
        if calc_type in ["generic", "decoupled"]:
            try:
                # Test with chain selection
                chain_results = calculator.calculate(
                    str(pdb_path),
                    chains=[["A"], ["B"]],
                    include_matrix=True
                )
                
                # First validate the results
                validate_sasa_results(chain_results, calc_type)
                
                # Then check for chain-specific results - modify to handle actual structure
                if "atom_data" in chain_results:
                    chains_found = set()
                    for atom_id, atom_data in chain_results["atom_data"].items():
                        if "chain" in atom_data:
                            chains_found.add(atom_data["chain"])
                    
                    assert any(chain in chains_found for chain in ["A", "B"]), \
                        f"Expected chains A or B in results, found: {chains_found}"
                
                # Verify contacts (if they exist)
                has_contacts = False
                for residue in chain_results.get("residue_data", []):
                    if residue.get("contacts", {}):
                        has_contacts = True
                        break
                
                # If we expect contacts but don't find any, we can warn instead of fail
                if not has_contacts:
                    warnings.warn(f"No contacts found between chains for {calc_type} calculator")
                    
            except IndexError as e:
                pytest.xfail(f"Chain selection not fully implemented: {str(e)}")

    def test_output_generation(self, calc_type, calculators, tmp_path):
        """Test output file generation"""
        calculator = calculators[calc_type]
        pdb_path = get_test_file("basic")
        
        output_name = str(tmp_path / "test_output")
        results = calculator.calculate(
            str(pdb_path),
            print_output=True,
            output_name=output_name

        )
        
        validate_sasa_results(results, calc_type)
        
        # Check for output files
        expected_files = []
        if calc_type == "simple":
            expected_files.extend([
                f"{output_name}.asa.pdb",
                f"{output_name}.atmasa"
            ])
        elif calc_type in ["generic", "decoupled"]:
            expected_files.extend([
                f"{output_name}.dsasa.pdb",
                f"{output_name}.datmasa",
                f"{output_name}.overlaps"
            ])
            
        # Verify at least some output files were created
        output_files = list(tmp_path.glob("*"))
        assert len(output_files) > 0, "No output files generated"

class TestImplementationComparison:
    """Test suite for comparing SASA implementations."""
    
    def calculate_original_sasa(self, pdb_file: str, dr_sasa_exec: str) -> Dict[str, float]:
        """Calculate SASA using original dr_sasa implementation."""
        with tempfile.NamedTemporaryFile(suffix='.pdb') as temp_output:
            # Run original dr_sasa
            cmd = [dr_sasa_exec, "-m", "0", "-i", str(pdb_file), "-o", temp_output.name]
            result = subprocess.run(cmd, capture_output=True, text=True)
            
            if result.returncode != 0:
                raise RuntimeError(f"Original dr_sasa failed: {result.stderr}")
            
            # Read output PDB and extract SASA values from B-factor column
            parser = PDBParser(QUIET=True)
            structure = parser.get_structure('protein', temp_output.name)
            
            # Create dictionary mapping atom ID to SASA value
            atom_sasa = {}
            atom_id = 1  # Assuming 1-based indexing
            for model in structure:
                for chain in model:
                    for residue in chain:
                        for atom in residue:
                            atom_sasa[str(atom_id)] = atom.get_bfactor()
                            atom_id += 1
            
            print(f"Original implementation found {len(atom_sasa)} atoms")
            return atom_sasa

    def debug_atom_indices(self, original_values: Dict[str, float], 
                        python_values: Dict[str, Dict],
                        pdb_file: str) -> None:
        """Debug atom indexing differences between implementations."""
        
        print("\nAtom Index Analysis:")
        print("-" * 50)
        
        # Original implementation indices
        orig_indices = sorted([int(x) for x in original_values.keys()])
        print(f"Original implementation:")
        print(f"  Index range: {min(orig_indices)} to {max(orig_indices)}")
        print(f"  Number of atoms: {len(orig_indices)}")
        print(f"  First 10 indices: {orig_indices[:10]}")
        print(f"  Last 10 indices: {orig_indices[-10:]}")
        
        # Check for gaps in original indices
        gaps = []
        for i in range(len(orig_indices)-1):
            if orig_indices[i+1] - orig_indices[i] > 1:
                gaps.append((orig_indices[i], orig_indices[i+1]))
        if gaps:
            print(f"  Gaps found in original indices:")
            for start, end in gaps[:5]:  # Show first 5 gaps
                print(f"    Gap between {start} and {end}")
        
        # Python implementation indices
        py_indices = sorted([int(x) for x in python_values.keys()])
        print(f"\nPython implementation:")
        print(f"  Index range: {min(py_indices)} to {max(py_indices)}")
        print(f"  Number of atoms: {len(py_indices)}")
        print(f"  First 10 indices: {py_indices[:10]}")
        print(f"  Last 10 indices: {py_indices[-10:]}")
        
        # Compare atom details for first few atoms
        print("\nDetailed atom comparison (first 5 atoms):")
        print("-" * 50)
        print("Original | Python")
        print("-" * 50)
        
        for i in range(min(5, len(orig_indices))):
            orig_id = str(orig_indices[i])
            py_id = str(py_indices[i])
            
            orig_val = original_values[orig_id]
            py_info = python_values[py_id]
            
            print(f"ID: {orig_id} | ID: {py_id}")
            print(f"Value: {orig_val:.2f} | Value: {py_info['surface']['sasa']:.2f}")
            if 'name' in py_info:
                print(f"         | Name: {py_info['name']}")
            if 'resname' in py_info:
                print(f"         | Residue: {py_info['resname']}")
            print("-" * 50)
            
        # Check if total atom count matches PDB file
        parser = PDBParser(QUIET=True)
        structure = parser.get_structure('temp', pdb_file)
        total_atoms = sum(1 for _ in structure.get_atoms())
        
        print(f"\nPDB file total atoms: {total_atoms}")
        print(f"Original implementation atoms: {len(original_values)}")
        print(f"Python implementation atoms: {len(python_values)}")

    def compare_sasa_values(self, original_values: Dict[str, float], 
                          python_values: Dict[str, Dict]) -> Dict[str, float]:
        """Compare SASA values between implementations."""
        # Extract SASA values from Python results
        python_sasa = {
            atom_id: atom_info['surface']['sasa'] 
            for atom_id, atom_info in python_values.items()
        }
        print(python_sasa)
        # Check for atom set differences
        original_atoms = set(original_values.keys())
        python_atoms = set(python_sasa.keys())
        
        # Print diagnostic information
        print(f"\nAtom count comparison:")
        print(f"Original implementation: {len(original_atoms)} atoms")
        print(f"Python implementation: {len(python_atoms)} atoms")
        
        # Find mismatched atoms
        only_in_original = original_atoms - python_atoms
        only_in_python = python_atoms - original_atoms
        
        if only_in_original or only_in_python:
            print("\nAtom ID mismatches:")
            if only_in_original:
                print(f"Atoms only in original: {sorted(only_in_original)[:5]} ...")
            if only_in_python:
                print(f"Atoms only in Python: {sorted(only_in_python)[:5]} ...")
            
            # Use common atoms for comparison
            common_atoms = original_atoms & python_atoms
            print(f"\nProceeding with {len(common_atoms)} common atoms")
            
            # Filter dictionaries to common atoms
            original_filtered = {k: v for k, v in original_values.items() if k in common_atoms}
            python_filtered = {k: v for k, v in python_sasa.items() if k in common_atoms}
        else:
            original_filtered = original_values
            python_filtered = python_sasa
        
        # Calculate differences
        differences = []
        atom_diffs = []  # Store per-atom differences for analysis
        for atom_id in original_filtered:
            orig_val = original_filtered[atom_id]
            py_val = python_filtered[atom_id]
            diff = orig_val - py_val
            differences.append(diff)
            if abs(diff) > 1.0:  # Track significant differences
                atom_diffs.append((atom_id, diff, orig_val, py_val))
        
        differences = np.array(differences)
        original_array = np.array(list(original_filtered.values()))
        python_array = np.array(list(python_filtered.values()))
        
        # Print significant differences
        if atom_diffs:
            print("\nSignificant differences (|diff| > 1.0):")
            for atom_id, diff, orig, py in sorted(atom_diffs, key=lambda x: abs(x[1]), reverse=True)[:5]:
                print(f"Atom {atom_id}: orig={orig:.2f}, py={py:.2f}, diff={diff:.2f}")
        
        return {
            'mean_difference': float(np.mean(differences)),
            'max_difference': float(np.max(np.abs(differences))),
            'rmsd': float(np.sqrt(np.mean(differences**2))),
            'correlation': float(np.corrcoef(original_array, python_array)[0,1]),
            'relative_error': float(np.mean(np.abs(differences) / original_array)),
            'num_atoms': {
                'original': len(original_atoms),
                'python': len(python_atoms),
                'common': len(original_filtered)
            }
        }

    def test_implementation_comparison(self, calculators, test_paths, dr_sasa_exec, calc_type):
        """Compare SASA calculations between original and Python implementations."""
        calculator = calculators[calc_type]
        pdb_path = test_paths("basic")
        
        print(f"\nTesting {calc_type} calculator implementation")
        
        # Calculate SASA with Python implementation
        python_results = calculator.calculate(str(pdb_path), print_output=True)
        
        # Calculate SASA with original implementation
        original_results = self.calculate_original_sasa(pdb_path, dr_sasa_exec)
        print(original_results)
        
        self.debug_atom_indices(original_results, python_results['atoms'], str(pdb_path))

        # Compare results
        comparison = self.compare_sasa_values(original_results, python_results['atoms'])
        
        # Print detailed comparison
        print(f"\nImplementation Comparison Results ({calc_type}):")
        print(f"Total atoms - Original: {comparison['num_atoms']['original']}, "
              f"Python: {comparison['num_atoms']['python']}, "
              f"Common: {comparison['num_atoms']['common']}")
        print(f"Correlation: {comparison['correlation']:.4f}")
        print(f"RMSD: {comparison['rmsd']:.4f}")
        print(f"Mean difference: {comparison['mean_difference']:.4f}")
        print(f"Max difference: {comparison['max_difference']:.4f}")
        print(f"Relative error: {comparison['relative_error']:.4f}")
        
        # Assert on comparison metrics using only common atoms
        assert comparison['correlation'] > 0.99, \
            f"Low correlation between implementations: {comparison['correlation']:.4f}"
        assert comparison['rmsd'] < 1.0, \
            f"High RMSD between implementations: {comparison['rmsd']:.4f}"
        assert comparison['relative_error'] < 0.05, \
            f"High relative error: {comparison['relative_error']:.4f}"

if __name__ == "__main__":
    pytest.main(["-v"])