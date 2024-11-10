# Python Bindings for dr_sasa_n (unofficial)

Python bindings for dr_sasa_n (Solvent Accessible Surface Area Calculator). If you use dr_sasa_n in your research, please cite:  
Ribeiro, J., Ríos-Vera, C., Melo, F. and Schüller, A. (2018) "Calculation of accurate contact surface areas between atoms for the quantitative analysis of non-bonded molecular interactions". Bioinformatics

## Overview

dr_sasa_n is a high-performance tool for calculating Solvent Accessible Surface Area (SASA) with the following features:

### Current Features
- High-performance SASA calculations on Linux systems with OpenMP support
- Efficient Python/C++ interface using pybind11
- NumPy-based memory management
- Multiple calculation modes for different analysis needs

### In Development
- Optional file I/O handling
- OpenCL/CUDA acceleration support
- Refined inter/intra BSA matrix calculations
- Comprehensive contact surface analysis

### Requirements
- Python 3.8+
- NumPy
- C++17 compatible compiler
- CMake 3.15+
- pybind11

### Important Notes



 **Asymmetric Contacts**: Contact surfaces are inherently asymmetric - CS(A→B) ≠ CS(B→A) 
 - Example: Large atom contacting small atom * A might lose 30Å² contacting B * B might only lose 15Å² contacting A 2. 
 
 **Result Files**: 
 - `*.asa.pdb`: PDB with SASA values in B-factor column 
 - `*.atmasa`: Detailed atom-by-atom analysis 
 - `*_vs_*.tsv`: Contact matrices between chains 
 - `*.overlaps`: Detailed overlap information

### Quick Install
```bash
# Clone repository
git clone https://github.com/username/dr_sasa_python.git
cd dr_sasa_python

# Install using provided script
./install.sh
```

## 1. SimpleSASA
Basic SASA calculator for single structures or complexes.

```python
from dr_sasa_py import SimpleSASA

# Initialize calculator
calculator = SimpleSASA(
    probe_radius=1.4,  # Water probe radius in Angstroms
    compute_mode=0     # CPU mode (GPU support coming soon)
)

# Calculate SASA
results = calculator.calculate(
    "protein.pdb",
    print_output=True,      # Generate detailed output files
    output_name="results"   # Base name for output files
)

# Access results
for atom_id, atom_data in results.items():
    print(f"Atom {atom_id}: SASA = {atom_data['sasa']:.2f} Å²")
```


## 2. GenericSASA
Advanced calculator for analyzing interactions between chain groups.

```python
from dr_sasa_py import GenericSASA

# Initialize calculator
calculator = GenericSASA(probe_radius=1.4)

# Calculate with chain selection
results = calculator.calculate(
    "complex.pdb",
    chains=[["A"], ["B"]],  # Analyze chains A and B
    include_matrix=True     # Generate interaction matrices
)

# Analyze chain interactions
for atom_id, atom_data in results.items():
    if atom_data['contacts']:
        print(f"\nAtom {atom_id} ({atom_data['chain']}) contacts:")
        for contact_id, contact_info in atom_data['contacts'].items():
            print(f"  - Contact with {contact_id}: {contact_info['contact_area']:.2f} Å²")
```

## 3. DecoupledSASA
Specialized calculator for complex molecular assemblies.

```python
from dr_sasa_py import DecoupledSASA

# Initialize calculator not working yet
calculator = DecoupledSASA(probe_radius=1.4)

# Calculate with multiple chain groups
results = calculator.calculate(
    "complex.pdb",
    chains=[["A", "B"], ["C", "D"]],  # Analyze AB vs CD interactions
    include_matrix=True
)
```

## Contributing
Contributions are welcome!
## License
This project is licensed under the MIT License - see the LICENSE file for details.
