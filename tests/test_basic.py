import pytest
import os, sys
from pathlib import Path

# Try different potential build paths
build_paths = [
    Path("build/lib"),
    Path("build/lib.linux-x86_64-cpython-310"),  # Common Linux path
    Path("build/lib.linux-x86_64-3.10"),         # Alternative Linux path
    Path("build/lib.macosx-10.9-x86_64-3.10"),   # Mac path
    Path("build"),                               # Direct build
]

# Add all potential paths
for build_path in build_paths:
    if build_path.exists():
        sys.path.append(str(build_path.absolute()))
        print(f"Added build path: {build_path.absolute()}")

def test_import():
    """Test different ways to import the module"""
    import_methods = []
    
    # Method 1: Direct import
    try:
        import dr_sasa_py
        import_methods.append("Direct import successful")
    except ImportError as e:
        print(f"Direct import failed: {e}")

    # Method 2: From import
    try:
        from dr_sasa_py import SimpleSASA
        import_methods.append("From import SimpleSASA successful")
    except ImportError as e:
        print(f"From import SimpleSASA failed: {e}")

    # Method 3: Import using importlib
    try:
        import importlib
        dr_sasa = importlib.import_module('dr_sasa_py')
        import_methods.append("Importlib import successful")
    except ImportError as e:
        print(f"Importlib import failed: {e}")

    # Method 4: Try importing with full path
    try:
        import os
        module_path = os.path.join(os.path.dirname(__file__), '..', 'build', 'lib')
        if module_path not in sys.path:
            sys.path.append(module_path)
        import dr_sasa_py as dr_sasa_full_path
        import_methods.append("Full path import successful")
    except ImportError as e:
        print(f"Full path import failed: {e}")

    # Print successful methods
    print("\nSuccessful import methods:")
    for method in import_methods:
        print(f"✓ {method}")

    # Assert at least one method worked
    assert len(import_methods) > 0, "No import methods succeeded"
    
    # Return the first successful import for use in other tests
    if "Direct import successful" in import_methods:
        import dr_sasa_py
        return dr_sasa_py
    elif "From import SimpleSASA successful" in import_methods:
        from dr_sasa_py import SimpleSASA
        return SimpleSASA
    elif "Importlib import successful" in import_methods:
        import importlib
        return importlib.import_module('dr_sasa_py')
    elif "Full path import successful" in import_methods:
        return dr_sasa_full_path

def test_basic_sasa():
    """Test basic SASA calculation with a PDB file"""
    # Get the module from our import test
    dr_sasa_module = test_import()
    
    # Create calculator using the successfully imported module
    if hasattr(dr_sasa_module, 'SimpleSASA'):
        calc = dr_sasa_module.SimpleSASA(probe_radius=1.4)
    else:
        calc = dr_sasa_module(probe_radius=1.4)
    
    # Calculate SASA
    test_pdb = os.path.join(os.path.dirname(__file__), "data", "3i40.pdb")
    assert os.path.exists(test_pdb), f"Test PDB file not found: {test_pdb}"
    
    results = calc.calculate_sasa(test_pdb)
    # Basic checks
    assert 'total_sasa' in results, "Missing total_sasa in results"
    assert results['total_sasa'] > 0, "Total SASA should be positive"
    assert 'atom_sasa' in results, "Missing atom_sasa in results"
    assert len(results['atom_sasa']) > 0, "No atoms found in results"
    
    print(f"\nTest Results:")
    print(f"Total SASA: {results['total_sasa']:.2f} Å²")
    print(f"Number of atoms: {len(results['atom_sasa'])}")
    print(f"First atom SASA: {results['atom_sasa'][0]:.2f} Å²")
    
    print("\nFirst 5 atoms:")
    for i in range(min(5, len(results['atom_sasa']))):
        print(f"Atom {results['atom_names'][i]} ({results['residue_names'][i]} {results['residue_numbers'][i]}): "
              f"{results['atom_sasa'][i]:.2f} Å²")
        
if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s"])  # -s flag to show print statements
