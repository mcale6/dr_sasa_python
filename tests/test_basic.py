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

TEST_DATA_DIR = os.path.join(os.path.dirname(__file__), "data")

import dr_sasa_py

def test_generic_sasa_basic():
    """Basic test of GenericSASA with default settings"""
    # Initialize calculator
    calculator = dr_sasa_py.GenericSASA(probe_radius=1.4)
    
    pdb_path = os.path.join(TEST_DATA_DIR, "3i40.pdb")
    if not os.path.exists(pdb_path):
        pytest.skip(f"Test PDB file not found: {pdb_path}")
        
    # Calculate with empty chains for automatic mode
    results = calculator.calculate(pdb_path, chains=[], include_matrix=False)
    
    # Basic validation
    assert isinstance(results, dict)
    assert 'atoms' in results
    assert 'surface' in results
    assert 'interactions' in results

def test_generic_sasa_advanced():
    calculator = dr_sasa_py.GenericSASA(probe_radius=1.4)
    
    pdb_path = os.path.join(TEST_DATA_DIR, "3i40.pdb")
    if not os.path.exists(pdb_path):
        pytest.skip(f"Test PDB file not found: {pdb_path}")
    
    # Use simple string lists for chains
    results = calculator.calculate(
        pdb_path,
        chains=[["A"], ["B"]],  # Simple nested list of strings
        include_matrix=False
    )
    
    # Alternative approach for automatic mode
    #auto_results = calculator.calculate(pdb_path, include_matrix=True)
    
    
    assert isinstance(results, dict)
    assert 'surface' in results
    assert np.sum(results['surface']['sasa']) > 0

if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s"])