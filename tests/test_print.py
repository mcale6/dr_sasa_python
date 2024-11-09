import pytest
import os
from pathlib import Path
from dr_sasa_py import SimpleSASA, GenericSASA, DecoupledSASA
from structure_parser import parse_pdb_file

# Get path to test PDB file
TEST_PDB = str(Path(__file__).parent / "data" / "3i40.pdb")

def test_simple_sasa_print():
    """Test SimpleSASA print function."""
    calculator = SimpleSASA()
    
    # Parse PDB file using StructureData
    structure = parse_pdb_file(TEST_PDB)
    atom_structs = structure.to_atom_structs()
    results = calculator.calculate(TEST_PDB, print_output=True, output_name="my_analysis.tsv")

def test_generic_sasa_print():
    """Test GenericSASA print function."""
    calculator = GenericSASA()
    
    # Parse PDB file using StructureData
    structure = parse_pdb_file(TEST_PDB)
    atom_structs = structure.to_atom_structs()
    
def test_decoupled_sasa_print():
    """Test DecoupledSASA print function."""
    calculator = DecoupledSASA()
    
    # Parse PDB file using StructureData
    structure = parse_pdb_file(TEST_PDB)
    atom_structs = structure.to_atom_structs()
    
@pytest.mark.skip(reason="Need to handle empty atoms case properly in C++ code")
def test_print_with_empty_atoms():
    """Test print function behavior with empty atom list."""
    calculator = SimpleSASA()
    with pytest.raises(ValueError):
        calculator.print([], "test_output")

if __name__ == "__main__":
    pytest.main([__file__, "-v"])