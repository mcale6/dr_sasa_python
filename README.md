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

### DecoupledSolver (mode 4)
- Contact points: $P_k = R_i * S_k + C_i - C_j$
- Pattern matching: $P(G) = \prod_{i \in G} p_i$ (primes)
- Buried area: $A_{buried} = A_{total} * (N_{pattern}/N_{total})$
- Contact area: $A_{contact} = A_{total} * (N_{contact}/N_{total})$
- Overlap area: $A_{overlap}(G) = A_{total} * (N_{overlap}/N_{total})$

The key difference is that GenericSolver computes overlaps by comparing states (original vs. new) while DecoupledSolver directly calculates overlaps within a single state, though both can output similar overlap and contact information in their results. DecoupledSolver doesnt not check wathever or not the atom is solvent exposed. 

## In Development
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

### `atoms`
Contains atom-level information with each atom ID as key:
```python
{
    "atom_id": {
        "name": str,            # Atom name
        "resname": str,         # Residue name
        "chain": str,           # Chain identifier
        "resid": int,          # Residue number
        "index": int,          # Atom index
        "coords": tuple,       # (x, y, z) coordinates
        
        "surface": {
            "sphere_area": float,    # Total sphere area
            "sasa": float,           # Solvent accessible surface area
            "buried_area": float,    # Buried surface area
            "contact_area": float,   # Contact surface area
            "dsasa": float,          # Delta SASA
            "total_asa": float,      # Total accessible surface area
            
            "backbone": {
                "total": float,      # Total backbone surface
                "polar": float,      # Polar backbone surface
                "hydrophobic": float # Hydrophobic backbone surface
            },
            "sidechain": {
                "total": float,      # Total sidechain surface
                "polar": float,      # Polar sidechain surface
                "hydrophobic": float # Hydrophobic sidechain surface
            },
            "groove": {              # For nucleic acids
                "major": {
                    "total": float,
                    "polar": float,
                    "hydrophobic": float
                },
                "minor": {
                    "total": float,
                    "polar": float,
                    "hydrophobic": float
                },
                "none": {
                    "total": float,
                    "polar": float,
                    "hydrophobic": float
                }
            },
            "ligand": {             # For ligands
                "total": float,
                "polar": float,
                "hydrophobic": float
            }
        },
        
        "properties": {
            "vdw": float,          # Van der Waals radius
            "polar": bool,         # Is polar
            "charge": str,         # Charge
            "struct_type": str     # Structure type
        },
        
        "contacts": {
            "nonbonded": {         # Non-bonded contacts
                "atom_id": {
                    "area": float,
                    "distance": float
                }
            },
            "overlap_groups": [     # Overlap information
                {
                    "atoms": list[int],
                    "area": float,
                    "normalized_area": float,
                    "buried_area": float
                }
            ]
        }
    }
}
```

### `residues`
Contains residue-level information with residue ID as key:
```python
{
    "residue_id": {
        "identifiers": {
            "chain": str,
            "name": str,
            "number": int
        },
        "structure": {
            "center": tuple,     # (x, y, z)
            "n_atoms": int
        },
        "surface": {
            "total_sasa": float,
            "total_area": float,
            "standard_sasa": float,
            "dsasa": float
        },
        "contacts": dict,        # Similar to atom contacts
        "overlaps": list        # Similar to atom overlaps
    }
}
```

### `chains`
Contains chain-level information with chain ID as key:
```python
{
    "chain_id": {
        "type": str,           # Chain type (PROTEIN, DNA, etc.)
        "residues": list[str]  # List of residue IDs
    }
}
```

### `lookup`
Contains cross-reference information:
```python
{
    "atoms": {
        "by_residue": dict,    # Residue ID to list of atom IDs
        "by_chain": dict       # Chain ID to list of atom IDs
    },
    "residues": {
        "by_chain": dict       # Chain ID to list of residue IDs
    }
}
```

### `metadata`
Contains analysis metadata:
```python
{
    "total_atoms": int,        # Total number of atoms
    "total_residues": int,     # Total number of residues
    "chains": list[str],       # List of chain IDs
    "types": list[str],        # List of molecule types
    "probe_radius": float      # Probe radius used
}
```

## Accessing the Data

### Atom Data
Access atom information with detailed surface and contact analysis:
```python
# Get data for atom with ID 1
atom = results["atoms"]["1"]

# Basic properties
print(f"Atom {atom['name']} in {atom['resname']}{atom['resid']}")
print(f"Location: x={atom['coords'][0]}, y={atom['coords'][1]}, z={atom['coords'][2]}")

# Surface analysis
surface = atom["surface"]
print(f"SASA: {surface['sasa']}")
print(f"Buried Area: {surface['buried_area']}")
print(f"Contact Area: {surface['contact_area']}")

# Component analysis
backbone = surface["backbone"]
print(f"Backbone surface: total={backbone['total']}, polar={backbone['polar']}")

# For nucleic acids, check groove surfaces
if atom['properties']['struct_type'] in ['DNA', 'RNA']:
    groove = surface["groove"]
    major = groove["major"]
    minor = groove["minor"]
    print(f"Major groove: total={major['total']}, polar={major['polar']}")
    print(f"Minor groove: total={minor['total']}, polar={minor['polar']}")
```

### Residue Data
Access residue-level information with surface and contact details:
```python
# Get residue by ID
res_id = "A_ALA_1"  # Format: "chain_resname_resid"
residue = results["residues"][res_id]

# Basic information
print(f"Residue: {residue['identifiers']['chain']}_{residue['identifiers']['name']}{residue['identifiers']['number']}")
print(f"Number of atoms: {residue['structure']['n_atoms']}")

# Surface analysis
surface = residue["surface"]
print(f"Total SASA: {surface['total_sasa']}")
print(f"Standard SASA: {surface['standard_sasa']}")
print(f"dSASA: {surface['dsasa']}")

# Contact analysis
contacts = residue["contacts"]
print(f"Number of contacts: {len(contacts)}")
for atom_id, contact in contacts.items():
    print(f"Contact with atom {atom_id}: area={contact['contact_area']}, distance={contact['distance']}")
```

### Chain Data
Access chain-level information and aggregated statistics:
```python
# Get chain information
chain = results["chains"]["A"]
print(f"Chain type: {chain['type']}")
print(f"Number of residues: {len(chain['residues'])}")

# List all residues in chain
for res_id in chain["residues"]:
    print(f"Residue: {res_id}")
```

### Using Lookup Tables
Efficiently navigate between atoms, residues, and chains:
```python
# Get all atoms in a residue
res_id = "A_ALA_1"
atom_ids = results["lookup"]["atoms"]["by_residue"][res_id]
for atom_id in atom_ids:
    atom = results["atoms"][atom_id]
    print(f"Atom {atom['name']}: SASA={atom['surface']['sasa']}")

# Get all residues in a chain
chain_id = "A"
res_ids = results["lookup"]["residues"]["by_chain"][chain_id]
for res_id in res_ids:
    residue = results["residues"][res_id]
    print(f"Residue {res_id}: {residue['surface']['total_sasa']}")
```

### Analyzing Surface Components
Access detailed surface component analysis:
```python
# Get backbone vs sidechain analysis for a protein atom
atom = results["atoms"]["1"]
if atom["properties"]["struct_type"] == "PROTEIN":
    backbone = atom["surface"]["backbone"]
    sidechain = atom["surface"]["sidechain"]
    print(f"Backbone: polar={backbone['polar']}, hydrophobic={backbone['hydrophobic']}")
    print(f"Sidechain: polar={sidechain['polar']}, hydrophobic={sidechain['hydrophobic']}")

# Get groove analysis for nucleic acid
if atom["properties"]["struct_type"] in ["DNA", "RNA"]:
    groove = atom["surface"]["groove"]
    for location in ["major", "minor", "none"]:
        g = groove[location]
        print(f"{location.capitalize()} groove: total={g['total']}, polar={g['polar']}, hydrophobic={g['hydrophobic']}")
```

### Metadata Access
Get analysis parameters and summary information:
```python
metadata = results["metadata"]
print(f"Total atoms analyzed: {metadata['total_atoms']}")
print(f"Total residues: {metadata['total_residues']}")
print(f"Chains present: {metadata['chains']}")
print(f"Molecule types: {metadata['types']}")
print(f"Probe radius used: {metadata['probe_radius']}")
```

## Quick Start Examples

### 1. Simple SASA Calculation
Basic SASA calculation for a structure:
```python
from dr_sasa import SimpleSASA

# Initialize calculator
calculator = SimpleSASA(probe_radius=1.4)
results = calculator.calculate("3i40.pdb")

# Access results
print(f"Total atoms analyzed: {results['metadata']['total_atoms']}")
```

### 2. Generic Analysis (Chain Interactions)
Analyze interactions between specific chains:
```python
from dr_sasa import GenericSASA

# Initialize calculator
calculator = GenericSASA(probe_radius=1.4)
results = calculator.calculate("3i40.pdb", chains=[["A"], ["B"]], print_output=True)

Selected complex surface (A^2):	3362.89
Object A complexed surface (A^2):	1201.52
Object B complexed surface (A^2):	2161.36
Object A uncomplexed surface (A^2):	1998.48
Object B uncomplexed surface (A^2):	2905.35
A <--- B buried surface (A^2):	796.957
A ---> B buried surface (A^2):	743.99
Interface A/B (A^2):	770.474
```

### 3. Decoupled Surface Analysis
Detailed surface component analysis:
```python
from dr_sasa import DecoupledSASA
from utils import convert_to_dataframes
# Decoupled does also check non-solvent exposed contacts
calc = dr_sasa_py.DecoupledSASA(probe_radius=1.4) 
result = calc.calculate(str("3i40.pdb"), chains=[["A"], ["B"]], # chains need to be defined!
    include_matrix=True, 
    print_output=True)
dfs = convert_to_dataframes(convert_to_dataframes)

# Access surface components
atoms = dfs['atoms']
print("\nSurface Composition:")
print(f"Backbone total: {atoms['backbone_total'].sum():.2f}")
print(f"Sidechain total: {atoms['sidechain_total'].sum():.2f}")
print("\nPolar/Hydrophobic Distribution:")
print(f"Polar surface: {atoms['polar_asa'].sum():.2f}")
print(f"Hydrophobic surface: {atoms['hyd_asa'].sum():.2f}")

A <--- B buried surface (A^2):	4787.38
A ---> B buried surface (A^2):	4910.12
Interface A/B (A^2):	4848.75
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

### Surface Point Calculation
- Points on atom i: p(k) = R(i) * s(k) + c(i)
  * where:
    - p(k) is the transformed point k
    - R(i) is radius of atom i
    - s(k) is unit sphere point k
    - c(i) is center of atom i

### Contact Check
- For each point against atom j: |p(k) - c(j)|² ≤ R(j)²
  * where:
    - p(k) is the point to check
    - c(j) is center of atom j
    - R(j) is radius of atom j

### 1. Point Generation
1. Start with unit sphere points
   - Generate points where sqrt(x² + y² + z²) = 1
2. Scale points by atom radius
   - point = point * R(I)
3. Translate to atom center
   - point = point + C(I)
4. For overlap checking, shift relative to other atom
   - point = point - C(J)

### 2. Burial Check
1. Calculate distance to other atom
   - dist(j) = (x - x(j))² + (y - y(j))² + (z - z(j))²
2. Check if point is buried
   - buried = true if dist(j) ≤ R(J)²

### 3. SASA Calculation
- Final SASA = (number of unburied points / total points) * 4πr²
  * where r is the atom radius plus probe radius

## Notes
- Points are typically distributed uniformly on a unit sphere
- The more points used, the more accurate the calculation
- Probe radius is added to atom radius for solvent accessibility


## Contributing
Contributions are welcome!

## License
This project is licensed under the MIT License - see the LICENSE file for details.
