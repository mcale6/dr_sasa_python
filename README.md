# DR_SASA Python Bindings

Python bindings for dr_sasa (Solvent Accessible Surface Area Calculator)

## Overview
dr_sasa is a high-performance tool for calculating Solvent Accessible Surface Area (SASA) and various forms of differential SASA analyses. These Python bindings provide efficient access to all DR_SASA functionality with minimal overhead.

## Installation

### Requirements
- Python 3.7+
- NumPy
- A C++17 compatible compiler
- CMake 3.15+

### Installation Steps



## Usage

### Available Modes

1. **Simple SASA (Mode 0)**
   - Basic SASA calculation for molecular structures
   ```python
   from drsasa import SimpleSASA
   
   calculator = SimpleSASA()
   calculator.set_probe_radius(1.4)  # Water probe radius in Angstroms
   result = calculator.run("protein.pdb")
   print(f"Total SASA: {result.statistics['total_sasa']}")
   ```

2. **Generic Delta SASA (Mode 1)**
   - Calculates SASA differences between chain groups
   ```python
   from drsasa import GenericDSASA
   
   calculator = GenericDSASA()
   calculator.set_chains([["A", "B"], ["C", "D"]])
   result = calculator.run("complex.pdb")
   ```

3. **Internal SASA Analysis**
   - **Residue Mode (Mode 2)**: Residue-level internal contacts
   - **Atom Mode (Mode 3)**: Atom-level internal contacts
   ```python
   from drsasa import InternalDSASA
   
   # Residue mode
   residue_calc = InternalDSASA.create_residue_mode()
   residue_result = residue_calc.run("protein.pdb")
   
   # Atom mode
   atom_calc = InternalDSASA.create_atom_mode()
   atom_result = atom_calc.run("protein.pdb")
   ```

4. **Decoupled Delta SASA (Mode 4)**
   - Separate surface calculations for molecular or chain contacts
   ```python
   from drsasa import DecoupledDSASA
   
   calculator = DecoupledDSASA()
   result = calculator.run("complex.pdb")
   ```

### Common Parameters

All calculators support these basic parameters:
```python
calculator.set_probe_radius(1.4)       # Set water probe radius
calculator.set_vdw_file("vdw.dat")     # Custom VdW radii
calculator.set_cl_mode(1)              # OpenCL acceleration mode
calculator.enable_reorder(True)        # Enable atom reordering
```

### Result Structure

All calculations return a `DRSASAResult` object containing:
- `sasa`: NumPy array of SASA values
- `rel_sasa`: Relative SASA values (when applicable)
- `coordinates`: Atomic coordinates
- `metadata`: Dictionary of additional data
- `statistics`: Calculation statistics

```python
result = calculator.run("protein.pdb")
print(result.sasa)           # SASA values
print(result.statistics)     # Statistical information
print(result.metadata)       # Additional metadata
```

## Advanced Features

### Matrix Output
Many modes support matrix output for detailed analysis:
```python
calculator.set_matrix_output(True)
result = calculator.run("protein.pdb")
matrix = result.statistics["interaction_matrix"]
```

### Chain Selection
Specify chain groups for analysis:
```python
calculator.set_chains([["A", "B"], ["C", "D"]])  # Compare AB vs CD
# or
calculator.set_chains([["A"]])  # Single chain analysis
```

### Quick Functions
Convenience functions for common operations:
```python
from drsasa import calculate_simple_sasa, calculate_delta_sasa

# Quick SASA calculation
result = calculate_simple_sasa("protein.pdb")

# Quick delta SASA between chains
result = calculate_delta_sasa("complex.pdb", [["A"], ["B"]])
```

## Performance Considerations

- Uses direct memory access with NumPy arrays
- Supports OpenCL acceleration
- Efficient handling of large structures
- Minimal Python/C++ boundary crossing
- No file writing and saving!

## Citation

If you use DR_SASA in your research, please cite:
```
Authors: Ribeiro J., Ríos-Vera C., Melo F., Schüller A.
Version: 0.5.0
```
