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
    assert "interaction_matrices" in results
    
    # Check intra matrices if present
    if "intra_matrices" in results:
        intra = results["intra_matrices"]
        assert all(k in intra for k in ['atom_matrix', 'residue_matrix', 'atom_labels'])
        assert isinstance(intra['atom_matrix'], np.ndarray)
        assert isinstance(intra['residue_matrix'], np.ndarray)

def test_decoupled_sasa_print():
    """Test DecoupledSASA print function and matrix generation."""
    calculator = DecoupledSASA()
    
    # Parse PDB file using StructureData
    structure = parse_pdb_file(TEST_PDB)
    atom_structs = structure.to_atom_structs()
    
    # Test with matrices
    results = calculator.calculate(TEST_PDB, print_output=True,
                                 include_matrix=True,
                                 output_name="decoupled_with_matrix.tsv")
    assert isinstance(results, dict)
    
    # Verify matrix structure if present
    if "interaction_matrices" in results:
        matrices = results["interaction_matrices"]
        assert isinstance(matrices, dict)
        
        # Check matrix content
        for key, data in matrices.items():
            if key != "residue_matrices":
                assert isinstance(data['matrix'], np.ndarray)
                assert data['matrix'].ndim == 2  # Should be 2D array

@pytest.mark.skip(reason="Need to handle empty atoms case properly in C++ code")
def test_print_with_empty_atoms():
    """Test print function behavior with empty atom list."""
    calculator = SimpleSASA()
    with pytest.raises(ValueError):
        calculator.print([], "test_output")

if __name__ == "__main__":
    pytest.main([__file__, "-v"])
