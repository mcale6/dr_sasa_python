import os
import sys
from pathlib import Path
build_path = Path("build")
sys.path.append(str(build_path.absolute()))
print(f"Added build path: {build_path.absolute()}")
import dr_sasa_py

import numpy as np

# Create and manipulate coordinates
coord = dr_sasa_py.Coord(1.0, 2.0, 3.0)
coord_array = coord.to_numpy()
new_coord = dr_sasa_py.Coord.from_numpy(np.array([4.0, 5.0, 6.0]))

# Create an atom
atom = dr_sasa_py.Atom()
atom.atom_name = "CA"
atom.res_name = "ALA"
atom.chain = "A"
atom.coords = coord
atom.radius = 1.7
atom.sasa = 0.0

# Create a residue
residue = dr_sasa_py.Residue()
residue.res_name = "ALA"
residue.chain = "A"
residue.res_index = 1
residue.atoms = [atom]

# Create and use VDW container
vdw = dr_sasa_py.VDWContainer("t.text")
vdw.generate_points()
points = vdw.get_points()  # Returns numpy array