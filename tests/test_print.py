import pytest
import os
import numpy as np
from pathlib import Path
from dr_sasa_py import SimpleSASA, GenericSASA, DecoupledSASA
from bindings.python.utils.structure_parser import StructureData, parse_pdb_file

# Get path to test PDB file
TEST_PDB = str(Path(__file__).parent / "data" / "3i40.pdb")

def test_generic_sasa_print():
    """Test GenericSASA print function and matrix generation."""
    calculator = GenericSASA()
    
    # Parse PDB file using StructureData
    structure = parse_pdb_file(TEST_PDB)
    atom_structs = structure.to_atom_structs()
    
    # Test with matrices
    results = calculator.calculate(TEST_PDB, print_output=True,
                                 include_matrix=True,
                                 output_name="generic_with_matrix.tsv")
    assert isinstance(results, dict)
    assert "inter_bsa_matrix" in results
    assert "intra_bsa_matrix" in results
    
    # Check intra matrices if present
    if "intra_matrices" in results:
        #assert all(k in intra for k in ['atom_matrix', 'residue_matrix', 'atom_labels'])
        assert isinstance(results["intra_bsa_matrix"]["atom_matrix"], np.ndarray)


@pytest.mark.skip(reason="Need to handle empty atoms case properly in C++ code")
def test_print_with_empty_atoms():
    """Test print function behavior with empty atom list."""
    calculator = SimpleSASA()
    with pytest.raises(ValueError):
        calculator.print([], "test_output")

if __name__ == "__main__":
    pytest.main([__file__, "-v"])
