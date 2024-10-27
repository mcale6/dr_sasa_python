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
    pdb_path = os.path.join(TEST_DATA_DIR, "PRED.pdb")
    results = calc.calculate(pdb_path)

    # Basic checks
    assert 'total_sasa' in results
    assert results['total_sasa'] > 0
    assert 'atom_sasa' in results
    assert len(results['atom_sasa']) > 0

    # Check consistency of values (not exact equality)
    total = sum(results['atom_sasa'])
    assert abs(total - results['total_sasa']) < 1e-6  # Allow small numerical differences

    print(f"\nSimpleSASA Results:")
    print(f"Total SASA: {results['total_sasa']:.2f}")
    print(f"Sum of atom SASA: {total:.2f}")
    print(f"Number of atoms: {len(results['atom_sasa'])}")

def test_decoupled_sasa():
    """Test DecoupledSASA class"""
    import dr_sasa_py

    # Test initialization
    calc = dr_sasa_py.DecoupledSASA(probe_radius=1.4, compute_mode=0)

    # Calculate SASA
    pdb_path = os.path.join(TEST_DATA_DIR, "PRED.pdb")
    results = calc.calculate(pdb_path)

    # Basic checks with detailed output
    print(f"\nDecoupledSASA Results:")
    print(f"Total SASA: {results['total_sasa']:.2f}")
    print(f"Number of atoms: {len(results['atom_sasa'])}")
    if results['total_sasa'] <= 0:
        print("WARNING: Zero or negative total SASA")
        print("First few atom SASAs:", results['atom_sasa'][:5])

    assert results['total_sasa'] > 0, "Total SASA should be positive"

def test_error_handling():
    """Test error handling in all classes"""
    import dr_sasa_py

    # Test with non-existent file
    with pytest.raises((RuntimeError, IOError, FileNotFoundError)):  # Accept multiple error types
        calc = dr_sasa_py.SimpleSASA()
        calc.calculate("nonexistent.pdb")

    # Test with invalid file content
    with pytest.raises((RuntimeError, ValueError)):
        # Create empty file
        with open("empty.pdb", "w") as f:
            f.write("")
        calc = dr_sasa_py.SimpleSASA()
        try:
            calc.calculate("empty.pdb")
        finally:
            os.remove("empty.pdb")  # Clean up

    # Test GenericSASA with invalid chains
    with pytest.raises((RuntimeError, ValueError)):
        calc = dr_sasa_py.GenericSASA()
        pdb_path = os.path.join(TEST_DATA_DIR, "PRED.pdb")
        calc.calculate(pdb_path, chains=[["X"]], mode=1)

def test_result_consistency():
    """Test consistency between different solvers"""
    import dr_sasa_py

    pdb_path = os.path.join(TEST_DATA_DIR, "PRED.pdb")

    # Calculate with all solvers
    simple = dr_sasa_py.SimpleSASA()
    results_simple = simple.calculate(pdb_path)

    generic = dr_sasa_py.GenericSASA()
    results_generic = generic.calculate(pdb_path, chains=[["A"]], mode=1)

    decoupled = dr_sasa_py.DecoupledSASA()
    results_decoupled = decoupled.calculate(pdb_path)

    print("\nConsistency Results:")
    print(f"SimpleSASA total: {results_simple['total_sasa']:.2f}")
    print(f"GenericSASA total: {results_generic['total_sasa']:.2f}")
    print(f"DecoupledSASA total: {results_decoupled['total_sasa']:.2f}")

    # Check for reasonable differences
    assert abs(results_simple['total_sasa'] - results_generic['total_sasa']) < results_simple['total_sasa'] * 0.1  # Allow 10% difference
    assert results_decoupled['total_sasa'] > 0

if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s"])  # -s flag to show print statements