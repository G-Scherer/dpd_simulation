# This file imports a GROMACS .top and .gro file and exports a lammps file

import MDAnalysis as mda
import numpy as np
import os

# Path to file
script_dir = os.path.dirname(os.path.abspath(__file__))

gro_file = os.path.join(script_dir, "final.gro")
top_file = os.path.join(script_dir, "setup.top")

# Build MDA universe
u = mda.Universe(top_file, gro_file, format="GRO", topology_format="ITP")

unique_types, type_ids = np.unique(u.atoms.types, return_inverse=True)

# Lammps expects stringified numbers and starts at 1
u.atoms.types = (type_ids + 1).astype(str)

# Write lammps file
lammps_out = os.path.join(script_dir, "lammps_raw.data")
with mda.Writer(lammps_out, format="DATA") as w:
    w.write(u.atoms)

print(f"LAMMPS file written sucessfully: {lammps_out}")