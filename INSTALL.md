# Installation Guide for dr_sasa_py

## Option 1: Using Conda (Recommended)

1. Create and activate conda environment:
```bash
# Create environment from yml file
conda env create -f environment.yml
conda activate dr_sasa

# Build the package
mkdir build
cd build
cmake ..
make
```

2. Add to PYTHONPATH:
```bash
export PYTHONPATH=$PYTHONPATH:$(pwd)/build/lib
```

## Option 2: Using Python venv

1. Create and activate virtual environment:
```bash
python -m venv .venv
source .venv/bin/activate  # On Windows: .venv\Scripts\activate
```

2. Install requirements:
```bash
pip install -r requirements.txt
```

3. Build the package:
```bash
mkdir build
cd build
cmake ..
make
```

4. Add to PYTHONPATH:
```bash
export PYTHONPATH=$PYTHONPATH:$(pwd)/build/lib
```

## Testing the Installation

```python
import dr_sasa_py

# Create calculator
calc = dr_sasa_py.SimpleSASA()

# Test with included PDB
results = calc.calculate_sasa("tests/data/3i40.pdb")
print(f"Total SASA: {results['total_sasa']:.2f}")
```

## Common Issues

1. If CMake can't find pybind11:
```bash
conda install -c conda-forge pybind11
# or
pip install pybind11
```

2. If the module isn't found after building:
- Make sure the build completed successfully
- Check that the .so file exists in build/lib
- Verify PYTHONPATH includes the build/lib directory