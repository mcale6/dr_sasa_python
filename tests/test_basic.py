import pytest
import os, sys
from pathlib import Path
build_path = Path("build")
sys.path.append(str(build_path.absolute()))
print(f"Added build path: {build_path.absolute()}")
from dr_sasa_py import DrSASA, ComputeBackend
import numpy as np
#python setup.py build_ext --inplace

# Define the path to test data
TEST_DATA_DIR = os.path.join(os.path.dirname(__file__), "data")

@pytest.fixture
def dr_sasa():
    """Create a DrSASA instance for testing"""
    return dr_sasa_py.DrSASA()

def test_initialization(dr_sasa):
    """Test basic initialization"""
    assert dr_sasa.probe_radius == 1.4
    assert dr_sasa.get_compute_device_info() is not None

def test_basic_sasa_calculation(dr_sasa):
    """Test basic SASA calculation"""
    pdb_path = os.path.join(TEST_DATA_DIR, "3i40.pdb")
    results = dr_sasa.calculate_sasa(pdb_path)
    
    # Basic checks
    assert 'total_sasa' in results
    assert results['total_sasa'] > 0
    assert 'atom_sasa' in results
    assert len(results['atom_sasa']) > 0
    
    # Check array consistency
    assert len(results['atom_sasa']) == len(results['atom_names'])
    assert len(results['coordinates']) == len(results['atom_names']) * 3

def test_delta_sasa_calculation(dr_sasa):
    """Test delta SASA calculation"""
    pdb_path = os.path.join(TEST_DATA_DIR, "3i40.pdb")
    results = dr_sasa.calculate_delta_sasa(
        pdb_path,
        chains=[["A"], ["B"]]
    )
    
    assert 'total_interface_area' in results
    assert results['total_interface_area'] > 0
    assert 'interface_atoms' in results

def test_parameter_setting(dr_sasa):
    """Test parameter setting"""
    dr_sasa.set_probe_radius(1.2)
    assert dr_sasa.probe_radius == 1.2
    
    dr_sasa.set_matrix_output(False)
    assert not dr_sasa.matrix_output

@pytest.mark.parametrize("probe_radius", [1.0, 1.4, 1.8])
def test_different_probe_radii(dr_sasa, probe_radius):
    """Test SASA calculation with different probe radii"""
    dr_sasa.set_probe_radius(probe_radius)
    pdb_path = os.path.join(TEST_DATA_DIR, "6gwp.pdb")
    results = dr_sasa.calculate_sasa(pdb_path)
    assert results['total_sasa'] > 0