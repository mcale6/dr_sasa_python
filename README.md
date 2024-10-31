# Python Bindings for dr_sasa_n

Python bindings for dr_sasa_n (Solvent Accessible Surface Area Calculator). If you use dr_sasa_n in your research, please cite: \
Ribeiro, J., Ríos-Vera, C., Melo, F. and Schüller, A. (2018) “Calculation of accurate contact surface areas between atoms for the quantitative analysis of non-bonded molecular interactions”. Bioinformatics

## Overview
dr_sasa_n is a high-performance tool for calculating Solvent Accessible Surface Area (SASA). +

- Works on linux system with openmp. 
- Minimal Python/C++ boundary crossing using pybinding11
- Uses memory access with NumPy arrays (maybe there is a better method)
In progress:
- file writing and saving as option
- Supports OpenCL/Cude acceleration 
- asa result table, check surface, overlaps results and add to utils. with atomidx.
- distribution analysis and inter/intra bsa matrix
- matrix is inconsistent, problem with indexing, might not be needed can be created separatly.

## Installation

### Requirements
- Python 3.8+
- NumPy
- A C++17 compatible compiler
- CMake 3.15+
- pybinding11

### Installation Steps

Check install.sh

# Available Modes

## Class Structure with pybindigs

### 1. SimpleSASA
Basic SASA calculator for single/complex structures.
compute_mode = 1 for gpu support, not supported yet.
```python
from dr_sasa_py import SimpleSASA

calculator = SimpleSASA(probe_radius=1.4, compute_mode=0) 
result = calculator.calculate("protein.pdb")
```

### 2. GenericSASA
Calculates delta SASA between chain groups with automatic mode detection. For protein chains <=2 use this.
```python
from dr_sasa_py import GenericSASA

calculator = GenericSASA(probe_radius=1.4, compute_mode=0)
# Modes automatically selected:
# (- Mode 4: Automatic (single chain group) not needed)
# - Mode 1: Manual (multiple chain groups)
# - Mode 5: Protein-protein (exactly two protein chains)
result = calculator.calculate("complex.pdb", chains=[["A"], ["B"]], include_matrix=True)
```

### 3. DecoupledSASA
Calculates contact surfaces between molecular components. When more chains involved are involved. (not tested)
```python
from dr_sasa_py import DecoupledSASA

calculator = DecoupledSASA(probe_radius=1.4, compute_mode=0)
# Modes automatically selected:
# - Mode 2: Molecular contacts
# - Mode 3: Chain contacts
result = calculator.calculate("complex.pdb", chains=[["AB", "CD"]])
```

## Notes

| Measurement | Mathematical Formula | Description |
|------------|---------------------|-------------|
| SASA | SASA(atom) = ∑(accessible points) × point_area | Solvent Accessible Surface Area for a single atom. Basic measure of exposed surface. |
| AREA_BURIED_BY_ATOM_area<br>(ov_table_area) | AREA_buried(atom_i, atom_j) = SASA(atom_i_alone) - SASA(atom_i_complex) | Raw buried area when atom i is in contact with atom j. Total area of overlap without normalization. |
| Normalized Area<br>(ov_norm_area) | AREA_norm(atom_i, overlap) = AREA_buried(atom_i) / N_atoms_in_overlap | Buried area normalized by number of atoms involved in the overlap. Distributes the buried area equally among all contributing atoms. |
| Contact Surface A->B | contact(A->B) = AREA_norm(A) | Surface area of atom A buried by contact with B. NOT equal to B->A contact! |
| Contact Surface B->A | contact(B->A) = AREA_norm(B) | Surface area of atom B buried by contact with A. |
| dSASA (Delta SASA) | dSASA = SASA(isolated) - SASA(complex) | Difference in accessible surface area between isolated and complexed states. |



 Normalized Area (ov_norm_area):
- Buried area divided by number of atoms in overlap
- Prevents double-counting in multi-atom contacts
- Gives fair distribution of burial contribution

Contact Surface:
- Bidirectional measure of atom interaction
- Sum of normalized contributions from both atoms
- Used for interface analysis
- Sum of CSA_ij atom = BSA

4. dSASA:
- Global measure of burial upon complexation
- Accounts for all changes in accessibility
- Used for binding interface analysis

Critical Point: Contact Surfaces are ASYMMETRIC!
- contact(A->B) ≠ contact(B->A)
- Example:
  * A large atom contacting a small atom
  * A might lose 30Å² of surface area due to contact with B
  * While B might only lose 15Å² of surface area due to contact with A.

This asymmetry is why the interaction matrices in the code are not symmetric:
```cpp
// These can be different values:
atom_matrix_AtoB[i][j]  // Surface of atom i buried by contact with atom j
atom_matrix_BtoA[j][i]  // Surface of atom j buried by contact with atom i
```

This is exactly why the code outputs two separate .tsv files for interactions:
1. A<-B: How much surface area atom A loses due to contact with B
2. B<-A: How much surface area atom B loses due to contact with A


## Result Structure

Check utils.py to see what results are returned.

## Example Analysis

```python

```