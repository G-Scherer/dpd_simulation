import MDAnalysis as mda
import numpy as np
import os

# Pfad zum aktuellen Verzeichnis
script_dir = os.path.dirname(os.path.abspath(__file__))

gro_file = os.path.join(script_dir, "final.gro")
top_file = os.path.join(script_dir, "setup.top")

u = mda.Universe(top_file, gro_file, format="GRO", topology_format="ITP")

unique_types, type_ids = np.unique(u.atoms.types, return_inverse=True)

u.atoms.types = (type_ids + 1).astype(str)

lammps_out = os.path.join(script_dir, "lammps1.data")
with mda.Writer(lammps_out, format="DATA") as w:
    w.write(u.atoms)

print(f"LAMMPS Datei erfolgreich erstellt: {lammps_out}")