# Python Bindings for dr_sasa_n (unofficial)

Python bindings for dr_sasa_n (Solvent Accessible Surface Area Calculator), a high-performance tool that calculates SASA using an extension of the Shrake-Ruply algorithm.

## Citation

If you use dr_sasa_n in your research, please cite:

> Ribeiro, J., Ríos-Vera, C., Melo, F. and Schüller, A. (2018) "Calculation of accurate contact surface areas between atoms for the quantitative analysis of non-bonded molecular interactions". Bioinformatics

Parameters for dSASA calculations are based on NACCESS (Chothia, 1976).

## Features

- High-performance SASA calculations with OpenMP/OpenCL support

- Efficient Python/C++ interface using pybind11

- NumPy-based memory management

- Multiple calculation modes:

  - Simple SASA

  - Generic SASA (chain-based analysis)

  - Decoupled SASA (independent contributions)

## Benchmark Dataset References

The implementation has been validated using datasets from:

1. Vangone, A. and Bonvin, A.M.J.J. (2015) "Contacts-based prediction of binding affinity in protein-protein complexes". eLife, e07454. DOI: 10.7554/eLife.07454

2. Moal, I.H., Agius, R., Bates, P.A. (2011) "Protein-protein binding affinity prediction on a diverse set of structures". Bioinformatics. DOI: 10.1093/bioinformatics/btr513


### Current Features
- High-performance SASA calculations on Linux systems with OpenMP/OpenCL support
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
When `print_output=True`, the following files are generated:

| File Extension | Description |
|----------------|-------------|
| `.asa.pdb` | PDB format file with SASA values stored in B-factor column |
| `.dasa.pdb` | PDB format file with delta-SASA values in B-factor column |
| `.atmasa` | Detailed per-atom SASA analysis |
| `.datmasa` | Detailed per-atom delta-SASA analysis |
| `_vs_.tsv` | Contact matrices between specified chains |
| `.overlaps` | Detailed atomic overlap information |
| `matrix.tsv` | Full atom-vs-atom contact matrices |

 **Result Dict**: 
The calculation returns a Python dictionary with the following structure:
When `include_matrix=True` and `print_output=True`, additional 3 fields are included ("residue_matrices", "interaction_matrices", "intra_matrices") beside atom_data and residue_data:

```python
{
    # Atom-level data
    "atom_data": {
        "1": {  # Atom ID as key
            "name": "CA",           # Atom name
            "resname": "ALA",       # Residue name
            "chain": "A",           # Chain identifier
            "resid": 1,            # Residue number
            "struct_type": "protein", # Structure type
            "coords": (23.456, 12.345, 34.567),  # XYZ coordinates
            "sphere_area": 120.5,   # Surface area
            "sasa": 45.6,          # Solvent accessible surface area
            "polar": True,         # Polarity flag
            "charge": 0.0          # Atomic charge
        },
        # ... more atoms
    },

    # Residue-level data
    "residue_data": [
        {
            "chain": "A",          # Chain identifier
            "resname": "ALA",      # Residue name
            "resid": 1,           # Residue number
            "total_sasa": 92.3,    # Total SASA for residue
            "total_area": 185.6,   # Total surface area
            "dsasa": 30.76,   # Difference from standard SASA
                                  # (STANDARD_SASA - total_sasa)
            "n_atoms": 5,          # Number of atoms in residue
            "center": (23.1, 12.8, 34.2),  # Center of mass
            
            # Contact information
            "contacts": {
                "2": {             # Contact atom ID
                    "struct_type": "protein",
                    "contact_area": 15.3,
                    "distance": 3.8
                },
                # ... more contacts
            },
            
            # Overlap information
            "overlaps": [
                {
                    "atoms": [2, 3, 4],     # Overlapping atom IDs
                    "overlap_area": 25.6,    # Area of overlap
                    "normalized_area": 0.212, # Normalized overlap area
                    "buried_area": 94.9      # Buried surface area
                },
                # ... more overlaps
            ]
        },
        # ... more residues
    ]
}
```

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
    compute_mode=0     # CPU mode (GPU support with OpenCL)
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

## 3. DecoupledSASA (not tested)
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
