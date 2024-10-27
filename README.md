# DR_SASA Python Bindings

Python bindings for the [dr_sasa_n](https://github.com/nioroso-x3/dr_sasa_n.git) library for calculating Solvent Accessible Surface Area.

! Still work in progress.

## Installation
*Clone Repo*
git clone --recursive https://github.com/mcale6/dr_sasa_python.git

*Dependencies*
pybind11 numpy pandas pytest cmake ocl-icd-opencl-dev

*Test import*
python dr_sasa_python/tests/test_import.py

Check install.sh for more details.

## Usage

Basic example:
```python
import dr_sasa_py

```

See the `examples/` directory for more usage examples.