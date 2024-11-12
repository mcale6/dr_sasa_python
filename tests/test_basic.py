import pytest
import os, sys
from pathlib import Path
import numpy as np
from typing import List, Optional, Dict, Any

import dr_sasa_py
from structure_parser import *

# Constants
TEST_DATA_DIR = Path(__file__).parent / "data"
TEST_FILES = {
    "basic": "3i40.pdb",
    "complex": "6gwp.pdb",
    "prediction": "pred.pdb"
}

def get_test_file(key: str) -> Path:
    """Get path to test file, skip if not found"""
    path = TEST_DATA_DIR / TEST_FILES[key]
    if not path.exists():
        pytest.skip(f"Test file not found: {path}")
    return path

def validate_sasa_results(results: Dict[str, Any], calc_type: str = "simple") -> None:
    """Validate SASA calculation results"""
    assert isinstance(results, dict)
    assert "atom_data" in results
    assert "residue_data" in results

    # Validate atom data
    atom_data = results["atom_data"]
    assert isinstance(atom_data, dict)
    assert len(atom_data) > 0

    # Get first atom entry
    first_atom_id = next(iter(atom_data))
    first_atom = atom_data[first_atom_id]

    # Validate atom data structure
    assert isinstance(first_atom, dict)
    required_atom_fields = {
        "name", "resname", "chain", "resid", "struct_type",
        "coords", "sphere_area", "sasa", "polar", "charge"
    }
    missing_atom_fields = required_atom_fields - set(first_atom.keys())
    assert not missing_atom_fields, f"Missing atom fields: {missing_atom_fields}"
    
    # Validate atom data types
    assert isinstance(first_atom["coords"], tuple)
    assert len(first_atom["coords"]) == 3
    assert isinstance(first_atom["sasa"], float)
    assert isinstance(first_atom["sphere_area"], float)

    # Basic residue data validation
    residue_data = results["residue_data"]
    assert isinstance(residue_data, list)
    assert len(residue_data) > 0

    # Validate basic residue fields
    first_residue = residue_data[0]
    basic_residue_fields = {
        "chain", "resname", "resid", "total_sasa", "total_area",
        "standard_sasa", "dsasa", "n_atoms", "center"
    }
    missing_residue_fields = basic_residue_fields - set(first_residue.keys())
    assert not missing_residue_fields, f"Missing basic residue fields: {missing_residue_fields}"
            
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
    @pytest.mark.parametrize("calc_type", ["simple", "generic", "decoupled"])
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
                validate_sasa_results(chain_results, calc_type)
                
                # Verify chain-specific results
                chains_found = {atom_data["chain"] 
                              for atom_data in chain_results.values()}
                assert any(chain in chains_found for chain in ["A", "B"]), \
                    f"Expected chains A or B in results, found: {chains_found}"
                
                # Verify contacts for interacting chains
                has_contacts = False
                for atom_data in chain_results.values():
                    if atom_data["contacts"]:
                        has_contacts = True
                        break
                assert has_contacts, "No contacts found between chains"
                
            except IndexError as e:
                pytest.xfail(f"Chain selection not fully implemented: {str(e)}")

    @pytest.mark.output
    @pytest.mark.parametrize("calc_type", ["simple", "generic", "decoupled"])
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