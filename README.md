# DR_SASA Python Bindings

Python bindings for the [dr_sasa_n](https://github.com/nioroso-x3/dr_sasa_n.git) library for calculating Solvent Accessible Surface Area.

! Still work in progress.

## Installation
Clone Repo: git clone --recursive https://github.com/mcale6/dr_sasa_python.git
Dependencies: pybind11 numpy pandas pytest cmake ocl-icd-opencl-dev
After installation test if the module can be imported
python dr_sasa_python/tests/test_import.py

Check install.sh for more details.

## Usage

Basic example:
```python
import dr_sasa_py

# Initialize
dr_sasa = dr_sasa_py.DrSASA()

# Calculate SASA
results = dr_sasa.calculate_sasa("example.pdb")
print(f"Total SASA: {results['total_sasa']:.2f}")
```

See the `examples/` directory for more usage examples.