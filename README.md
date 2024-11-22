# Python Bindings for dr_sasa_n (unofficial)

Python bindings for dr_sasa_n (Solvent Accessible Surface Area Calculator), a high-performance tool that calculates SASA using an extension of the Shrake-Ruply algorithm.
\
If you use dr_sasa_n in your research, please cite:

> Ribeiro, J., Ríos-Vera, C., Melo, F. and Schüller, A. (2018) "Calculation of accurate contact surface areas between atoms for the quantitative analysis of non-bonded molecular interactions". Bioinformatics

Parameters for dSASA calculations are based on NACCESS (Chothia, 1976).

## Features

- High-performance SASA calculations with OpenMP/OpenCL support

- Efficient Python/C++ interface using pybind11

- Optional original print output

- Returns Extensive Analysis Results Data Structure

- Multiple calculation modes: SimpleSolver > DecoupledSolver > GenericSolver (Computational Efficiency)


### SimpleSASA (mode 0)
- Basic SASA: $A_{SASA} = A_{total} * (N_{accessible}/N_{total})$

### GenericSolver (modes 1-3)
- $dSASA = SASA_{chainAB} - SASA_{chainA}$
- Pattern matching: $P(G) = \prod_{i \in G} p_i$ (primes)
- Buried area: $A_{buried} = A_{total} * (N_{pattern}/N_{total})$

### DecoupledSolver (mode 4)
- Contact points: $P_k = R_i * S_k + C_i - C_j$
- Contact area: $A_{contact} = A_{total} * (N_{contact}/N_{total})$
- Overlap area: $A_{overlap}(G) = A_{total} * (N_{overlap}/N_{total})$

The key difference is that GenericSolver computes overlaps by comparing states (original vs. new) while DecoupledSolver directly calculates overlaps within a single state, though both can output similar overlap and contact information in their results. DecoupledSolver doesnt not check wathever or not the atom is solvent exposed. 

### Point Transformations
$$
\begin{align*}
\text{Surface Points:}& \quad \vec{p}_k = R_i\vec{s}_k + \vec{c}_i \\
\text{Contact Check:}& \quad ||\vec{p}_k - \vec{c}_j||^2 \leq R_j^2
\end{align*}
$$

## In Development
- Fixing Decoupled SASA
- Input validation (e.g empty atoms)
- CUDA support
- Comprehensive contact surface analysis, plots

## Benchmark Dataset References

The implementation has been validated using datasets from:

1. Vangone, A. and Bonvin, A.M.J.J. (2015) "Contacts-based prediction of binding affinity in protein-protein complexes". eLife, e07454. DOI: 10.7554/eLife.07454

2. Moal, I.H., Agius, R., Bates, P.A. (2011) "Protein-protein binding affinity prediction on a diverse set of structures". Bioinformatics. DOI: 10.1093/bioinformatics/btr513

Check out [Benchmark Results](data/README.md)

## Quick Installation

### Requirements
- Python 3.8+
- NumPy
- C++17 compatible compiler
- CMake 3.15+
- pybind11

### One-Line Installation
Download and run the installation script:
```bash
curl -s https://raw.githubusercontent.com/mcale6/dr_sasa_python/main/install.sh | bash
```

### Manual Installation Steps

1. Install system dependencies:
```bash
sudo apt-get update
sudo apt-get install -y build-essential cmake git python3 python3-dev python3-venv python3-full ocl-icd-opencl-dev
```

2. Create and activate a virtual environment:
```bash
python3 -m venv ~/dr_sasa_venv
source ~/dr_sasa_venv/bin/activate
```

3. Install Python dependencies:
```bash
pip install --upgrade pip setuptools wheel
pip install "pybind11[global]" numpy pandas pytest
```

4. Clone and build the repository:
```bash
git clone --recursive https://github.com/mcale6/dr_sasa_python.git
cd dr_sasa_python
mkdir -p build && cd build
cmake ..
make -j4
```

5. Set up Python path:
```bash
echo "export PYTHONPATH=$(pwd)/build/lib:$(pwd)/build:$(pwd)" >> ~/.bashrc
source ~/.bashrc
```

### Verifying Installation
Test the installation by importing the package:
```bash
python -c "import dr_sasa_py; print('DR-SASA Python installed successfully!')"
```

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

# Analysis Results Data Structure
The analysis results are returned as a dictionary containing atom-level data, residue-level data, and a residue index for efficient lookups.\
When `include_matrix=True` and `print_output=True`, additional 3 fields are included beside atom_data and residue_data: "residue_matrices", "interaction_matrices", "intra_matrices". 

## Structure Overview
```python
{
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

    "residue_data": [  # List of residues in order
        {
            "chain": "A",          # Chain identifier
            "resname": "ALA",      # Residue name
            "resid": 1,           # Residue number
            "total_sasa": 92.3,    # Total SASA for residue
            "total_area": 185.6,   # Total surface area
            "dsasa": 30.76,        # Difference from standard SASA
            "n_atoms": 5,          # Number of atoms in residue
            "center": (23.1, 12.8, 34.2),  # Center of mass
            
            "contacts": {          # Contacts with other atoms
                "2": {             # Contact atom ID
                    "struct_type": "protein",
                    "contact_area": 15.3,
                    "distance": 3.8
                }
                # ... more contacts
            },
            
            "overlaps": [         # List of overlap groups
                {
                    "atoms": [2, 3, 4],     # Overlapping atom IDs
                    "overlap_area": 25.6,    # Area of overlap
                    "normalized_area": 0.212, # Normalized overlap area
                    "buried_area": 94.9      # Buried surface area
                }
                # ... more overlaps
            ]
        }
        # ... more residues
    ],

    "residue_index": {  # Maps residue identifiers to list indices
        "A_ALA_1": 0,   # Format: "chain_resname_resid": index
        "A_GLY_2": 1,
        # ... more residue indices
    }
}
```

## Accessing the Data

### Atom Data
Access atom information directly by atom ID:
```python
# Get data for atom with ID 1
atom = results["atom_data"]["1"]
print(f"Atom {atom['name']} has SASA of {atom['sasa']}")
```

### Residue Data
Two ways to access residue data:

1. Sequential access (by index):
```python
# Get first residue in structure
first_res = results["residue_data"][0]
print(f"First residue is {first_res['resname']}{first_res['resid']}")

# Iterate through all residues
for residue in results["residue_data"]:
    print(f"Residue {residue['chain']}_{residue['resname']}{residue['resid']} has {residue['n_atoms']} atoms")
```

2. Lookup by residue identifier:
```python
# Get specific residue using residue_index
res_id = "A_ALA_1"  # Format: "chain_resname_resid"
idx = results["residue_index"][res_id]
residue = results["residue_data"][idx]
print(f"SASA for {res_id}: {residue['total_sasa']}")
```

### Working with Contacts and Overlaps
```python
# Access contacts for a residue
residue = results["residue_data"][0]
for atom_id, contact in residue["contacts"].items():
    print(f"Contact with atom {atom_id}: area = {contact['contact_area']}, distance = {contact['distance']}")

# Access overlaps for a residue
for overlap in residue["overlaps"]:
    print(f"Overlap involving atoms {overlap['atoms']}: area = {overlap['overlap_area']}")
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

## Sphere Point Distribution Comparison
Point distribution: $dA = r^2 \sin(\theta) d\theta d\phi$
### Golden Spiral
- Deterministic pattern
- Points follow spiral from pole to pole
- Equal area between points but not equal distances
- *Point density: ∝ 1/sin(θ), error: E_golden ∝ ∫(1/sin(θ)) dθ, variable with orientation*
- Error varies with burial orientation due to uneven distribution
- Higher accuracy near poles
- Lower accuracy near equator
- Not rotationally invariant

### Thomson
Thomson optimization: $E = \sum_{i}\sum_{j\neq i} \frac{1}{|\vec{r}_i - \vec{r}_j|}$
- Uniform error distribution
- Equal distances between neighboring points
- Rotationally invariant results
- More accurate overall surface area calculation
- *Point density: ≈ constant, error ≈ constant, independent of orientation*
- Better representation of buried surface patches

### SASA Calculation Process

$$
\begin{align*}
\text{1. Point Generation:} \\
&\bullet \text{ Unit sphere points: } \sqrt{x^2+y^2+z^2}=1 \\
&\bullet \text{ Scale: } \vec{p} = p \cdot R_I \\
&\bullet \text{ Translate: } \vec{p} = \vec{p} + \vec{C_I} \\
&\bullet \text{ For overlap: } \vec{p} = \vec{p} - \vec{C_J} \\
\\
\text{2. Burial Check:} \\
&\text{dist}_j = (x-x_j)^2 + (y-y_j)^2 + (z-z_j)^2 \\
&\text{buried} = (\text{dist}_j \leq R_J^2) \\
\\
\text{3. SASA Calculation:} \\
&SASA = \frac{N_{\text{unburied}}}{N_{\text{total}}} \cdot 4\pi r^2
\end{align*}
$$


## Contributing
Contributions are welcome!

## License
This project is licensed under the MIT License - see the LICENSE file for details.
