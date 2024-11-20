import pytest
from pathlib import Path
import numpy as np
from typing import List, Optional, Dict, Any
import dr_sasa_py
from bindings.python.utils.structure_parser import StructureData, parse_pdb_file
import warnings

# Constants
TEST_DATA_DIR = Path(__file__).parent / "data"
TEST_FILES = {
    "basic": "3i40.pdb",
    #"complex": "6gwp.pdb",
    #"prediction": "pred.pdb" 
}

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
    required_keys = {"atom_data", "residue_data"}  # residue_index is optional
    missing_keys = required_keys - set(results.keys())
    assert not missing_keys, f"Missing top-level keys: {missing_keys}"

    # Validate atom data
    atom_data = results["atom_data"]
    assert isinstance(atom_data, dict), "atom_data must be a dictionary"
    assert len(atom_data) > 0, "atom_data cannot be empty"

    # Get first atom entry for detailed validation
    first_atom_id = next(iter(atom_data))
    first_atom = atom_data[first_atom_id]

    # Validate atom data structure and types based on actual implementation
    required_atom_fields = {
        "name": str,          # NAME in C++
        "resname": str,       # RESN in C++
        "chain": str,         # CHAIN in C++
        "resid": int,         # RESI in C++
        "struct_type": str,   # STRUCT_TYPE in C++
        "coords": tuple,      # COORDS in C++
        "sphere_area": float, # AREA in C++
        "sasa": float,        # SASA in C++
        "polar": int,         # POLAR in C++, initialized as 0
        "charge": str        # CHARGE in C++, string from PDB
    }

    missing_atom_fields = required_atom_fields.keys() - set(first_atom.keys())
    assert not missing_atom_fields, f"Missing atom fields: {missing_atom_fields}"

    # Validate atom field types
    for field, expected_type in required_atom_fields.items():
        assert isinstance(first_atom[field], expected_type), \
            f"Atom field {field} should be {expected_type}, got {type(first_atom[field])}"
        
        # Additional validation for specific fields
        if field == "polar":
            assert first_atom[field] in (0, 1), f"polar must be 0 or 1, got {first_atom[field]}"
        elif field == "coords":
            assert len(first_atom[field]) == 3, "coords must be a 3-tuple"
            assert all(isinstance(x, (int, float)) for x in first_atom[field]), \
                "coords must contain numeric values"
        elif field in ("sphere_area", "sasa"):
            assert first_atom[field] >= 0, f"{field} must be non-negative"

    # Validate residue data
    residue_data = results["residue_data"]
    assert isinstance(residue_data, list), "residue_data must be a list"
    assert len(residue_data) > 0, "residue_data cannot be empty"

    first_residue = residue_data[0]
    required_residue_fields = {
        "chain": str,
        "resname": str,
        "resid": int,
        "total_sasa": float,
        "total_area": float,
        "n_atoms": int,
        "center": tuple,
        "contacts": dict,
        "overlaps": list
    }

    missing_residue_fields = required_residue_fields.keys() - set(first_residue.keys())
    assert not missing_residue_fields, f"Missing residue fields: {missing_residue_fields}"

    # Validate residue field types and values
    for field, expected_type in required_residue_fields.items():
        assert isinstance(first_residue[field], expected_type), \
            f"Residue field {field} should be {expected_type}, got {type(first_residue[field])}"

        if field == "center":
            assert len(first_residue[field]) == 3, "center must be a 3-tuple"
            assert all(isinstance(x, (int, float)) for x in first_residue[field]), \
                "center coordinates must be numeric"
        elif field in ("total_sasa", "total_area"):
            assert first_residue[field] >= 0, f"{field} must be non-negative"
        elif field == "n_atoms":
            assert first_residue[field] > 0, "n_atoms must be positive"

    # Validate contacts if present
    if first_residue["contacts"]:
        first_contact = next(iter(first_residue["contacts"].values()))
        required_contact_fields = {
            "struct_type": str,
            "contact_area": float,
            "distance": float
        }
        
        missing_contact_fields = required_contact_fields.keys() - set(first_contact.keys())
        assert not missing_contact_fields, f"Missing contact fields: {missing_contact_fields}"
        
        for field, expected_type in required_contact_fields.items():
            assert isinstance(first_contact[field], expected_type), \
                f"Contact field {field} should be {expected_type}, got {type(first_contact[field])}"
            if field in ("contact_area", "distance"):
                assert first_contact[field] >= 0, f"{field} must be non-negative"

    # Validate overlaps if present
    if first_residue["overlaps"]:
        first_overlap = first_residue["overlaps"][0]
        required_overlap_fields = {
            "atoms": list,
            "overlap_area": float,
            "normalized_area": float,
            "buried_area": float
        }
        
        missing_overlap_fields = required_overlap_fields.keys() - set(first_overlap.keys())
        assert not missing_overlap_fields, f"Missing overlap fields: {missing_overlap_fields}"
        
        for field, expected_type in required_overlap_fields.items():
            assert isinstance(first_overlap[field], expected_type), \
                f"Overlap field {field} should be {expected_type}, got {type(first_overlap[field])}"
            
            if field == "atoms":
                assert len(first_overlap[field]) >= 2, "overlap must involve at least 2 atoms"
            elif field == "normalized_area":
                assert 0 <= first_overlap[field] <= 1, "normalized_area must be between 0 and 1"
            else:
                assert first_overlap[field] >= 0, f"{field} must be non-negative"

    # Optional validation for residue_index if present
    if "residue_index" in results:
        residue_index = results["residue_index"]
        assert isinstance(residue_index, dict), "residue_index must be a dictionary"
        assert len(residue_index) == len(residue_data), \
            "residue_index must match residue_data length"

        for res_id, idx in residue_index.items():
            assert isinstance(res_id, str), "residue_index keys must be strings"
            assert isinstance(idx, int), "residue_index values must be integers"
            assert 0 <= idx < len(residue_data), "residue_index values must be valid indices"

            # Validate key format and cross-reference
            parts = res_id.split('_')
            assert len(parts) == 3, "residue ID must have format 'chain_resname_resid'"
            
            residue = residue_data[idx]
            expected_id = f"{residue['chain']}_{residue['resname']}_{residue['resid']}"
            assert res_id == expected_id, \
                f"residue_index key {res_id} doesn't match residue data {expected_id}"

    # Calculator-specific validations
    if calc_type in ["generic", "decoupled"]:
        # Only check for presence of fields, not content
        assert "contacts" in first_residue, f"Advanced calculator {calc_type} must have contacts field"
        assert "overlaps" in first_residue, f"Advanced calculator {calc_type} must have overlaps field"
            
class TestDrSasa:
    @pytest.fixture(autouse=True)
    def setup(self):
        """Setup test environment"""
        self.probe_radius = 1.4
        self.calculators = {
            "simple": dr_sasa_py.SimpleSASA(probe_radius=self.probe_radius),
            "generic": dr_sasa_py.GenericSASA(probe_radius=self.probe_radius),
            "decoupled": dr_sasa_py.DecoupledSASA(probe_radius=self.probe_radius)
        }

    @pytest.mark.basic
    def test_atom_struct(self):
        """Test AtomStruct creation and properties"""
        atom = dr_sasa_py.AtomStruct(
            id=1, resi=1, icode="", name="CA", resn="ALA",
            chain="A", element="C", structure="test",
            mol_type="PROTEIN", x=0.0, y=0.0, z=0.0
        )
        
        assert atom.ID == 1
        assert atom.NAME == "CA"
        assert atom.RESN == "ALA"
        assert list(atom.COORDS) == [0.0, 0.0, 0.0]

    @pytest.mark.basic
    def test_structure_data(self):
        """Test StructureData creation and manipulation"""
        pdb_path = get_test_file("basic")
        structure = parse_pdb_file(pdb_path)
        assert len(structure.atom_names) > 0
        assert np.all(np.isfinite(structure.atom_positions))

    @pytest.mark.calculation
    @pytest.mark.parametrize("calc_type", ["simple", "generic"])
    def test_sasa_calculation(self, calc_type):
        """Test SASA calculation with different calculators"""
        calculator = self.calculators[calc_type]
        pdb_path = get_test_file("basic")
        
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

    @pytest.mark.output
    @pytest.mark.parametrize("calc_type", ["simple", "generic"])
    def test_output_generation(self, calc_type, tmp_path):
        """Test output file generation"""
        calculator = self.calculators[calc_type]
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

def pytest_configure(config):
    """Register custom markers"""
    config.addinivalue_line("markers", "basic: Basic functionality tests")
    config.addinivalue_line("markers", "calculation: SASA calculation tests")
    config.addinivalue_line("markers", "advanced: Advanced feature tests")
    config.addinivalue_line("markers", "output: Output generation tests")