import pytest
import dr_sasa_py
import numpy as np

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
    results = dr_sasa.calculate_sasa("tests/data/1bl0.pdb")
    
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
    results = dr_sasa.calculate_delta_sasa(
        "tests/data/3i40.pdb",
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
    results = dr_sasa.calculate_sasa("tests/data/6gwp.pdb")
    assert results['total_sasa'] > 0