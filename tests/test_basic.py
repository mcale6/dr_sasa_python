import pytest
import os, sys
from pathlib import Path
import numpy as np

# Try different potential build paths
build_paths = [
    Path("build/lib"),
    Path("build/lib.linux-x86_64-cpython-310"),  
    Path("build/lib.linux-x86_64-3.10"),         
    Path("build/lib.macosx-10.9-x86_64-3.10"),   
    Path("build"),                               
]

# Add all potential paths
for build_path in build_paths:
    if build_path.exists():
        sys.path.append(str(build_path.absolute()))
        print(f"Added build path: {build_path.absolute()}")

import dr_sasa_py

def test_basic_sasa():
    """Test basic SASA calculation with and without matrices"""
    calc = dr_sasa_py.SimpleSASA(probe_radius=1.4)
    test_pdb = os.path.join(os.path.dirname(__file__), "data", "3i40.pdb")
    assert os.path.exists(test_pdb), f"Test PDB file not found: {test_pdb}"
    
    # Test without matrices
    results_no_matrix = calc.calculate(test_pdb, include_matrix=False)
    assert 'matrices' not in results_no_matrix, "Matrices present when disabled"
    
    # Test with matrices
    results = calc.calculate(test_pdb, include_matrix=True)
    
    # Basic checks
    assert 'atom_info' in results, "Missing atom_info in results"
    assert 'sasa' in results, "Missing sasa in results"
    
    # Matrix checks if present
    if 'matrices' in results:
        matrices = results['matrices']
        assert 'inter_molecular' in matrices, "Missing inter-molecular matrices"
        assert 'intra_molecular' in matrices, "Missing intra-molecular matrices"
        
        # Print matrix info
        print("\nMatrix information:")
        print(f"Inter-molecular matrices present: {list(matrices['inter_molecular'].keys())}")
        print(f"Intra-molecular matrices present: {list(matrices['intra_molecular'].keys())}")

if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s"])