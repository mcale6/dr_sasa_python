import pytest
import dr_sasa_py

def test_residue_analysis(dr_sasa):
    """Test residue-level analysis"""
    results = dr_sasa.calculate_residue_sasa(
        "tests/data/3i40.pdb",
        chain=["A"]
    )
    
    assert 'residues' in results
    assert len(results['residues']) > 0
    
    # Check residue data structure
    first_residue = results['residues'][0]
    assert 'chain' in first_residue
    assert 'number' in first_residue
    assert 'name' in first_residue
    assert 'sasa' in first_residue

def test_contact_surface(dr_sasa):
    """Test contact surface calculation"""
    results = dr_sasa.calculate_contact_surface(
        "tests/data/3i40.pdb",
        chain=["A", "B"]
    )
    
    assert 'contacts' in results
    assert 'total_contact_area' in results
    assert results['total_contact_area'] > 0

@pytest.mark.parametrize("backend", [
    dr_sasa_py.ComputeBackend.CPU,
    pytest.param(dr_sasa_py.ComputeBackend.OPENCL, 
                marks=pytest.mark.skipif(not hasattr(dr_sasa_py, 'OPENCL'),
                                      reason="OpenCL not available"))
])
def test_compute_backends(backend):
    """Test different compute backends"""
    dr_sasa = dr_sasa_py.DrSASA(backend=backend)
    results = dr_sasa.calculate_sasa("tests/data/3i40.pdb")
    assert results['total_sasa'] > 0