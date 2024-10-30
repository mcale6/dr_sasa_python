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
    """Test basic SASA calculation with detailed result checking"""
    calc = dr_sasa_py.SimpleSASA(probe_radius=1.4)
    test_pdb = os.path.join(os.path.dirname(__file__), "data", "3i40.pdb")
    assert os.path.exists(test_pdb), f"Test PDB file not found: {test_pdb}"
    
    # Test without matrices first
    results_no_matrix = calc.calculate(test_pdb, include_matrix=False)
    assert 'matrices' not in results_no_matrix, "Matrices present when disabled"
    
    # Test with matrices
    results = calc.calculate(test_pdb, include_matrix=True)
    
    # Basic checks (existing)
    assert 'atom_info' in results, "Missing atom_info in results"
    assert 'sasa' in results, "Missing sasa in results"
    assert 'molecular' in results, "Missing molecular in results"
    
    # Matrix checks
    assert 'matrices' in results, "Missing matrices in results"
    matrices = results['matrices']
    
    # Check inter-molecular matrices
    assert 'inter_molecular' in matrices, "Missing inter-molecular matrices"
    inter = matrices['inter_molecular']
    required_inter_fields = [
        'atomic', 'residue', 
        'col_atoms', 'row_atoms',
        'col_types', 'row_types'
    ]
    for field in required_inter_fields:
        assert field in inter, f"Missing {field} in inter-molecular matrices"
    
    # Check intra-molecular matrices
    assert 'intra_molecular' in matrices, "Missing intra-molecular matrices"
    intra = matrices['intra_molecular']
    required_intra_fields = [
        'atomic', 'residue',
        'col_atoms', 'row_atoms',
        'col_types', 'row_types'
    ]
    for field in required_intra_fields:
        assert field in intra, f"Missing {field} in intra-molecular matrices"
    
    # Verify matrix properties
    def check_matrix_properties(matrix_data, name):
        assert len(matrix_data['atomic']) > 0, f"Empty atomic matrix in {name}"
        assert len(matrix_data['residue']) > 0, f"Empty residue matrix in {name}"
        assert len(matrix_data['col_atoms']) == len(matrix_data['row_atoms']), \
            f"Mismatched dimensions in {name} matrix"
    
    check_matrix_properties(inter, "inter-molecular")
    check_matrix_properties(intra, "intra-molecular")
    
    # Print matrix information
    print("\nMatrix Analysis:")
    print("Inter-molecular matrices:")
    print(f"  Atomic matrix size: {len(inter['atomic'])}")
    print(f"  Residue matrix size: {len(inter['residue'])}")
    print(f"  Number of column atoms: {len(inter['col_atoms'])}")
    print(f"  Number of row atoms: {len(inter['row_atoms'])}")
    
    print("\nIntra-molecular matrices:")
    print(f"  Atomic matrix size: {len(intra['atomic'])}")
    print(f"  Residue matrix size: {len(intra['residue'])}")
    print(f"  Number of column atoms: {len(intra['col_atoms'])}")
    print(f"  Number of row atoms: {len(intra['row_atoms'])}")

def test_matrix_calculations():
    """Test specific matrix generation and properties"""
    calc = dr_sasa_py.GenericSASA(probe_radius=1.4)
    test_pdb = os.path.join(os.path.dirname(__file__), "data", "3i40.pdb")
    
    # Test with multiple chains
    results = calc.calculate(test_pdb, [["A"], ["B"]], include_matrix=True)
    matrices = results['matrices']
    
    # Check matrix consistency
    inter = matrices['inter_molecular']
    assert all(isinstance(val, float) for val in inter['atomic']), \
        "Non-float values in atomic matrix"
    assert all(isinstance(val, float) for val in inter['residue']), \
        "Non-float values in residue matrix"
    
    # Check matrix values
    assert all(val >= 0 for val in inter['atomic']), \
        "Negative values in atomic matrix"
    assert all(val >= 0 for val in inter['residue']), \
        "Negative values in residue matrix"
    
    # Test matrix dimensions
    n_atoms = len(results['atom_info']['names'])
    n_residues = len(set(zip(results['atom_info']['mol_types'], 
                            [i//3 for i in range(n_atoms*3)])))  # Approximate
    
    print("\nMatrix Dimensions:")
    print(f"Number of atoms: {n_atoms}")
    print(f"Estimated number of residues: {n_residues}")
    print(f"Atomic matrix entries: {len(inter['atomic'])}")
    print(f"Residue matrix entries: {len(inter['residue'])}")

def test_all_calculators():
    """Test matrix generation in all calculator types"""
    test_pdb = os.path.join(os.path.dirname(__file__), "data", "3i40.pdb")
    
    calculators = [
        (dr_sasa_py.SimpleSASA(), {}),
        (dr_sasa_py.GenericSASA(), {"chains": [["A"], ["B"]]}),
        (dr_sasa_py.DecoupledSASA(), {"chains": [["A"], ["B"]]}),
        (dr_sasa_py.RelativeSASA(), {})
    ]
    
    for calc, kwargs in calculators:
        # Test with matrices
        results = calc.calculate(test_pdb, include_matrix=True, **kwargs)
        assert 'matrices' in results, f"Missing matrices in {calc.__class__.__name__}"
        
        # Test without matrices
        results_no_matrix = calc.calculate(test_pdb, include_matrix=False, **kwargs)
        assert 'matrices' not in results_no_matrix, \
            f"Matrices present when disabled in {calc.__class__.__name__}"

if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s"])