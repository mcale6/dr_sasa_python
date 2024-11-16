import pytest
import os, sys
from pathlib import Path

def test_import():
    """Test different ways to import the module"""
    import_methods = []
    
    # Method 1: Direct import
    try:
        import dr_sasa_py
        import_methods.append("Direct import successful")
    except ImportError as e:
        print(f"Direct import failed: {e}")

    # Method 2: From import all SASA classes
    try:
        from dr_sasa_py import SimpleSASA, GenericSASA, DecoupledSASA
        import_methods.append("From import SASA classes successful")
    except ImportError as e:
        print(f"From import SASA classes failed: {e}")

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

if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s"])  # -s flag to show print statements