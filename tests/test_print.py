import pytest
from dr_sasa_py import SimpleSASA, GenericSASA, DecoupledSASA
import numpy as np

@pytest.fixture
def sample_atoms():
    # Create a calculator and get atoms from a small test PDB
    calculator = SimpleSASA()
    return calculator.calculate("tests/data/3i40.pdb")

def test_simple_sasa_print(sample_atoms):
    calculator = SimpleSASA()
    output = calculator.print(sample_atoms, "simple_test")
    assert isinstance(output, str)
    assert len(output) > 0
    # Verify output format
    assert "REMARK 000 SASA SOLVER" in output

def test_generic_sasa_print(sample_atoms):
    calculator = GenericSASA()
    output = calculator.print(sample_atoms, "generic_test")
    assert isinstance(output, str)
    assert len(output) > 0
    # Check format specific to generic SASA output
    assert any("ATOM" in line or "HETATM" in line for line in output.split('\n'))

def test_decoupled_sasa_print(sample_atoms):
    calculator = DecoupledSASA()
    output = calculator.print(sample_atoms, "decoupled_test")
    assert isinstance(output, str)
    assert len(output) > 0
    # Check format specific to decoupled SASA output
    assert any("ATOM" in line or "HETATM" in line for line in output.split('\n'))

def test_print_with_empty_atoms():
    calculator = SimpleSASA()
    empty_atoms = []  # or however you represent empty atoms
    with pytest.raises(Exception):  # adjust exception type as needed
        calculator.print(empty_atoms, "empty_test")

def test_print_with_invalid_fname():
    calculator = SimpleSASA()
    with pytest.raises(Exception):  # adjust exception type as needed
        calculator.print(sample_atoms, "")  # or other invalid filename cases

if __name__ == "__main__":
    pytest.main([__file__])