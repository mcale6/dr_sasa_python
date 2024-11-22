import pytest
import os
import numpy as np
from pathlib import Path
from dr_sasa_py import SimpleSASA, GenericSASA, DecoupledSASA
from bindings.python.utils.structure_parser import StructureData, parse_pdb_file

# Get path to test PDB file
TEST_PDB = str(Path(__file__).parent / "data" / "3i40.pdb")


@pytest.mark.skip(reason="Need to handle empty atoms case properly in C++ code")
def test_print_with_empty_atoms():
    """Test print function behavior with empty atom list."""
    calculator = SimpleSASA()
    with pytest.raises(ValueError):
        calculator.print([], "test_output")

if __name__ == "__main__":
    pytest.main([__file__, "-v"])
