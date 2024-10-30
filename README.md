# dr_sasa Python Bindings

Python bindings for dr_sasa_n (Solvent Accessible Surface Area Calculator). If you use dr_sasa_n in your research, please cite: \
Ribeiro, J., Ríos-Vera, C., Melo, F. and Schüller, A. (2018) “Calculation of accurate contact surface areas between atoms for the quantitative analysis of non-bonded molecular interactions”. Bioinformatics

## Overview
dr_sasa is a high-performance tool for calculating Solvent Accessible Surface Area (SASA). +

- Works on linux system with openmp. 
- Minimal Python/C++ boundary crossing using pybinding11
- Uses memory access with NumPy arrays (maybe there is a better method)
- In progress: file writing and saving as option

- In progress: Supports OpenCL/Cude acceleration 
- In progress: asa, interaction, overlaps and bsa table 

## Installation

### Requirements
- Python 3.7+
- NumPy
- A C++17 compatible compiler
- CMake 3.15+
- pybinding11

### Installation Steps

Check install.sh \
Under tests/run_dr_sasa_n.sh one can se example usage for the c++ code. (needs to be compiled on its own!)

## Usage

### Available Modes
# dr_sasa_n Python Bindings

Python bindings for dr_sasa (Solvent Accessible Surface Area Calculator)

## Class Structure with pybindigs

### 1. SimpleSASA
Basic SASA calculator for single structures.
compute_mode = 1 for gpu support, not supported yet.
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
result = calculator.calculate("complex.pdb", chains=[["A"], ["B"]], include_matrix=false)
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

```

## Example Analysis

```python

```