import pytest
import os
import sys
from pathlib import Path

# Add build path to Python path
build_path = Path("build/lib")
sys.path.append(str(build_path.absolute()))
print(f"Added build path: {build_path.absolute()}")

# Test data path
TEST_DATA_DIR = os.path.join(os.path.dirname(__file__), "data")

def test_simple_sasa():
    """Test SimpleSASA class"""
    import dr_sasa_py
    
    # Test initialization
    calc = dr_sasa_py.SimpleSASA(probe_radius=1.4, compute_mode=0)
    
    # Calculate SASA for test PDB
    pdb_path = os.path.join(TEST_DATA_DIR, "3i40.pdb")
    results = calc.calculate(pdb_path)
    
    # Basic checks
    assert 'total_sasa' in results
    assert results['total_sasa'] > 0
    assert 'atom_sasa' in results
    assert len(results['atom_sasa']) > 0
    
    # Check returned data structures
    assert len(results['atom_names']) == len(results['atom_sasa'])
    assert len(results['coordinates']) == len(results['atom_sasa']) * 3
    
    print(f"\nSimpleSASA Results:")
    print(f"Total SASA: {results['total_sasa']:.2f}")
    print(f"Number of atoms: {len(results['atom_sasa'])}")

def test_generic_sasa():
    """Test GenericSASA class with chain analysis"""
    import dr_sasa_py
    
    # Test initialization
    calc = dr_sasa_py.GenericSASA(probe_radius=1.4, compute_mode=0)
    
    # Calculate SASA with chain analysis
    pdb_path = os.path.join(TEST_DATA_DIR, "3i40.pdb")
    results = calc.calculate(pdb_path, chains=[["A"], ["B"]], mode=1)
    
    # Basic checks
    assert 'total_sasa' in results
    assert results['total_sasa'] > 0
    
    # Check chain-specific data
    assert 'sasa_by_type' in results
    
    print(f"\nGenericSASA Results:")
    print(f"Total SASA: {results['total_sasa']:.2f}")
    if 'sasa_by_type' in results:
        for mol_type, sasa in results['sasa_by_type'].items():
            print(f"SASA for {mol_type}: {sasa:.2f}")

def test_decoupled_sasa():
    """Test DecoupledSASA class"""
    import dr_sasa_py
    
    # Test initialization
    calc = dr_sasa_py.DecoupledSASA(probe_radius=1.4, compute_mode=0)
    
    # Calculate SASA
    pdb_path = os.path.join(TEST_DATA_DIR, "3i40.pdb")
    results = calc.calculate(pdb_path)
    
    # Basic checks
    assert 'total_sasa' in results
    assert results['total_sasa'] > 0
    
    print(f"\nDecoupledSASA Results:")
    print(f"Total SASA: {results['total_sasa']:.2f}")

def test_error_handling():
    """Test error handling in all classes"""
    import dr_sasa_py
    
    # Test with non-existent file
    with pytest.raises(Exception):
        calc = dr_sasa_py.SimpleSASA()
        calc.calculate("nonexistent.pdb")
    
    # Test with invalid chains
    with pytest.raises(Exception):
        calc = dr_sasa_py.GenericSASA()
        calc.calculate(os.path.join(TEST_DATA_DIR, "3i40.pdb"), 
                      chains=[["X"]], mode=1)  # Non-existent chain

@pytest.mark.parametrize("probe_radius", [1.0, 1.4, 1.8])
def test_probe_radius(probe_radius):
    """Test different probe radii with all classes"""
    import dr_sasa_py
    
    pdb_path = os.path.join(TEST_DATA_DIR, "3i40.pdb")
    
    # Test SimpleSASA
    simple = dr_sasa_py.SimpleSASA(probe_radius=probe_radius)
    results_simple = simple.calculate(pdb_path)
    
    # Test GenericSASA
    generic = dr_sasa_py.GenericSASA(probe_radius=probe_radius)
    results_generic = generic.calculate(pdb_path, chains=[["A"]], mode=1)
    
    # Test DecoupledSASA
    decoupled = dr_sasa_py.DecoupledSASA(probe_radius=probe_radius)
    results_decoupled = decoupled.calculate(pdb_path)
    
    print(f"\nResults for probe radius {probe_radius}:")
    print(f"SimpleSASA: {results_simple['total_sasa']:.2f}")
    print(f"GenericSASA: {results_generic['total_sasa']:.2f}")
    print(f"DecoupledSASA: {results_decoupled['total_sasa']:.2f}")

def test_result_consistency():
    """Test consistency of results between different solvers"""
    import dr_sasa_py
    
    pdb_path = os.path.join(TEST_DATA_DIR, "3i40.pdb")
    
    # Calculate with all solvers
    simple = dr_sasa_py.SimpleSASA()
    results_simple = simple.calculate(pdb_path)
    
    generic = dr_sasa_py.GenericSASA()
    results_generic = generic.calculate(pdb_path, chains=[["A"]], mode=1)
    
    decoupled = dr_sasa_py.DecoupledSASA()
    results_decoupled = decoupled.calculate(pdb_path)
    
    # Compare number of atoms
    n_atoms = len(results_simple['atom_sasa'])
    assert len(results_generic['atom_sasa']) == n_atoms
    assert len(results_decoupled['atom_sasa']) == n_atoms
    
    print("\nConsistency check:")
    print(f"Number of atoms: {n_atoms}")
    print(f"SimpleSASA total: {results_simple['total_sasa']:.2f}")
    print(f"GenericSASA total: {results_generic['total_sasa']:.2f}")
    print(f"DecoupledSASA total: {results_decoupled['total_sasa']:.2f}")

if __name__ == "__main__":
    pytest.main([__file__, "-v"])