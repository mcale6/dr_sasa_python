# DR_SASA Python Bindings

Python bindings for dr_sasa (Solvent Accessible Surface Area Calculator). If you use DR_SASA in your research, please cite: \
Ribeiro, J., Ríos-Vera, C., Melo, F. and Schüller, A. (2018) “Calculation of accurate contact surface areas between atoms for the quantitative analysis of non-bonded molecular interactions”. Bioinformatics

## Overview
dr_sasa is a high-performance tool for calculating Solvent Accessible Surface Area (SASA). +

- Work on linux system with openmp.
- Minimal Python/C++ boundary crossing using pybinding11
- Uses memory access with NumPy arrays (maybe there is a better method)

- In progress: No file writing and saving optional
- In progress: Supports OpenCL/Cude acceleration 
- In progress: Contact analysis


## Installation

### Requirements
- Python 3.7+
- NumPy
- A C++17 compatible compiler
- CMake 3.15+
- pybinding11

### Installation Steps

Check install.sh


## Usage

### Available Modes
# DR_SASA Python Bindings

Python bindings for dr_sasa (Solvent Accessible Surface Area Calculator)

[Previous general sections remain the same until Usage...]

## Class Structure

### 1. SimpleSASA
Basic SASA calculator for single structures.
```python
from dr_sasa_py import SimpleSASA

calculator = SimpleSASA(probe_radius=1.4, compute_mode=0)
result = calculator.calculate("protein.pdb")
```

### 2. GenericSASA
Calculates delta SASA between chain groups with automatic mode detection.
```python
from dr_sasa_py import GenericSASA

calculator = GenericSASA(probe_radius=1.4, compute_mode=0)
# Modes automatically selected:
# - Mode 4: Automatic (single chain group)
# - Mode 1: Manual (multiple chain groups)
# - Mode 5: Protein-protein (exactly two protein chains)
result = calculator.calculate("complex.pdb", chains=[["A"], ["B"]])
```

### 3. DecoupledSASA
Calculates contact surfaces between molecular components.
```python
from dr_sasa_py import DecoupledSASA

calculator = DecoupledSASA(probe_radius=1.4, compute_mode=0)
# Modes automatically selected:
# - Mode 2: Molecular contacts
# - Mode 3: Chain contacts
result = calculator.calculate("complex.pdb", chains=[["A", "B"]])
```

### 4. RelativeSASA
Calculates relative accessibility compared to reference state.
```python
from dr_sasa_py import RelativeSASA

calculator = RelativeSASA(probe_radius=1.4, compute_mode=0)
result = calculator.calculate("protein.pdb")
```

## Result Structure

All calculations return a dictionary with the following structure:

### 1. Atom Information
```python
result["atom_info"] = {
    "names": list,          # Atom names (e.g., "CA", "N", "O")
    "elements": list,       # Element types
    "mol_types": list,      # Molecule types (e.g., "PROTEIN", "DNA")
    "coordinates": array,   # Nx3 array of atomic coordinates
    "radii": array,        # VdW radii
    "occupancies": array,  # Occupancy values
    "b_factors": array,    # Temperature factors
    "charges": list,       # Atomic charges
    "is_hetatm": list,    # HETATM flags
    "atom_types": array    # Internal atom type classifications
}
```

### 2. SASA Results
```python
result["sasa"] = {
    "values": array,       # Per-atom SASA values
    "delta": array,        # Delta SASA values (when applicable)
    "relative": array,     # Relative SASA values
    "by_mol_type": dict   # SASA summed by molecule type
}
```

### 3. Interface Analysis
```python
result["interface"] = {
    "total_area": float,              # Total interface area
    "buried_surface_area": float,     # Total buried surface
    "by_chain": dict,                 # Interface area per chain
    "atoms": dict                     # Interface atoms by chain
}
```

### 4. Molecular Analysis
```python
result["molecular"] = {
    "atoms_by_type": dict,     # Count of atoms per molecular type
    "residues_by_type": dict,  # Count of residues per type
    "type_contacts": dict      # Contacts between molecular types
}
```

### 5. Interaction Data (when available)
```python
result["interactions"] = {
    "by_type": dict,          # Interactions grouped by molecular types
    "energies": array         # Interaction energies
}
```

## Example Analysis

```python
calculator = SimpleSASA(probe_radius=1.4)
result = calculator.calculate("protein.pdb")

# Basic SASA analysis
total_sasa = result["sasa"]["values"].sum()
per_residue_sasa = {}
for i, resname in enumerate(result["atom_info"]["residue_names"]):
    per_residue_sasa[resname] = result["sasa"]["values"][i]

# Molecular composition
mol_composition = result["molecular"]["atoms_by_type"]
print(f"Molecule contains: {mol_composition}")

# Interface analysis (for applicable modes)
if "interface" in result:
    print(f"Total interface area: {result['interface']['total_area']}")
    print(f"Buried surface area: {result['interface']['buried_surface_area']}")
```