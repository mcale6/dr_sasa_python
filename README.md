# DR_SASA Python Bindings

Python bindings for the [dr_sasa_n](https://github.com/nioroso-x3/dr_sasa_n.git) library for calculating Solvent Accessible Surface Area.

Still work in progress.

## Installation

### Prerequisites
- Conda (or Miniconda)
- CMake (>=3.18)
- C++ compiler with C++17 support

### Installation Steps

1. Clone the repository:
```bash
git clone https://github.com/mcale6/dr_sasa_python.git
cd dr_sasa_python
```

2. Create and activate conda environment:
```bash
conda env create -f environment.yml
conda activate dr_sasa_env
```

3. Install the package:
```bash
pip install -e .
```

### Running Tests
```bash
pytest tests/
```

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