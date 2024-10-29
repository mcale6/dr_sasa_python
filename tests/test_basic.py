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
    # Get the module from our import test    
    # Create calculator
    if hasattr(dr_sasa_py, 'SimpleSASA'):
        calc = dr_sasa_py.SimpleSASA(probe_radius=1.4)
    else:
        calc = dr_sasa_py(probe_radius=1.4)
    
    # Calculate SASA
    test_pdb = os.path.join(os.path.dirname(__file__), "data", "3i40.pdb")
    assert os.path.exists(test_pdb), f"Test PDB file not found: {test_pdb}"
    
    results = calc.calculate(test_pdb)
    
    # 1. Check atom_info dictionary
    assert 'atom_info' in results, "Missing atom_info in results"
    atom_info = results['atom_info']
    required_atom_fields = [
        'names', 'elements', 'mol_types', 'coordinates', 
        'radii', 'occupancies', 'b_factors', 'charges', 
        'is_hetatm', 'atom_types'
    ]
    for field in required_atom_fields:
        assert field in atom_info, f"Missing {field} in atom_info"
        assert len(atom_info[field]) > 0, f"Empty {field} in atom_info"
    
    # 2. Check SASA results
    assert 'sasa' in results, "Missing sasa in results"
    sasa_info = results['sasa']
    required_sasa_fields = ['values', 'delta', 'relative', 'by_mol_type']
    for field in required_sasa_fields:
        assert field in sasa_info, f"Missing {field} in sasa_info"
    
    # 3. Check molecular information
    assert 'molecular' in results, "Missing molecular in results"
    mol_info = results['molecular']
    assert 'atoms_by_type' in mol_info, "Missing atoms_by_type"
    assert 'residues_by_type' in mol_info, "Missing residues_by_type"
    assert 'type_contacts' in mol_info, "Missing type_contacts"
    
    # 4. Check coordinates
    coords = atom_info['coordinates']
    assert len(coords.shape) == 2, "Coordinates should be 2D array"
    assert coords.shape[1] == 3, "Coordinates should have 3 columns (x,y,z)"
    
    # Print summary statistics
    print(f"\nTest Results Summary:")
    print(f"Number of atoms: {len(sasa_info['values'])}")
    print(f"Total SASA: {sasa_info['values'].sum():.2f} Å²")
    print(f"\nMolecule composition:")
    for mol_type, count in mol_info['atoms_by_type'].items():
        print(f"  {mol_type}: {count} atoms")
    
    print("\nSample atom details:")
    for i in range(min(5, len(atom_info['names']))):
        print(f"Atom {i+1}: {atom_info['names'][i]} ({atom_info['elements'][i]}) - "
              f"SASA: {sasa_info['values'][i]:.2f} Å²")
    
    # Optional: Check interface information if available
    if 'interface' in results:
        interface_info = results['interface']
        print("\nInterface information:")
        print(f"Total interface area: {interface_info['total_area']:.2f} Å²")
        print(f"Buried surface area: {interface_info['buried_surface_area']:.2f} Å²")
        
    # Optional: Check interaction information if available
    if 'interactions' in results:
        interaction_info = results['interactions']
        print("\nInteraction information available:")
        print(f"Number of interaction types: {len(interaction_info['by_type'])}")

def test_advanced_checks():
    """Additional tests for specific modes or calculations"""    
    # Test that values are reasonable
    calc = dr_sasa_py.SimpleSASA(probe_radius=1.4)
    test_pdb = os.path.join(os.path.dirname(__file__), "data", "3i40.pdb")
    results = calc.calculate(test_pdb)
    
    # Check value ranges
    sasa_values = results['sasa']['values']
    assert all(v >= 0 for v in sasa_values), "Found negative SASA values"
    assert all(v < 1000 for v in sasa_values), "Unreasonably large SASA values"
    
    # Check coordinate validity
    coords = results['atom_info']['coordinates']
    assert not np.any(np.isnan(coords)), "Found NaN in coordinates"
    assert not np.any(np.isinf(coords)), "Found Inf in coordinates"

if __name__ == "__main__":
    pytest.main([__file__, "-v", "-s"])