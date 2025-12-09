import argparse
import os
import numpy as np
import random
from scipy.constants import Boltzmann,Avogadro

parser = argparse.ArgumentParser(description='Create water system')

# define arguments, box default taken from past sim
parser.add_argument("--box", type=float, default=10, help="box length in nm")
parser.add_argument("--density", type=float, default=3, help="number density")
# 128: wpl:24, apl:0.85, 256: wpl: 22.38, apl: 0.88
#----------------
# parse arguments
args = parser.parse_args()

# save variables
box = args.box
density = args.density

r_ref = 0.711 #nm
e_ref = Boltzmann*298.15 #K
q = 8.861242189860825 

box = box/r_ref
box_x = box_y = box_z = box
box_vol = box**3
n_water = int(box_vol*density)

random.seed(75)

print(f"Creating system:")
print(f"  Water: {n_water}")
print(f"  Box volume: {box_vol} reduced units")
print(f" density: {n_water/box_vol}")

script_dir = os.path.dirname(os.path.abspath(__file__))

lines = []

# Add header and basic information about system

lines.append("Water setup")

lines.append("")

lines.append(f"{n_water} atoms")
lines.append(f"0 bonds")
lines.append(f"0 angles")
lines.append("0 dihedrals")
lines.append("0 impropers")

lines.append("")

lines.append("1 atom types")
lines.append("0 bond types")
lines.append("0 angle types")
lines.append("0 dihedral types")
lines.append("0 improper types")

lines.append("")

lines.append(f"0 {box_x} xlo xhi")
lines.append(f"0 {box_y} ylo yhi")
lines.append(f"0 {box_z} zlo zhi")

lines.append("")

lines.append("Masses")

lines.append("")

lines.append("1 1")

lines.append("")

lines.append("Atoms")

lines.append("")

for id in range(1,n_water+1):

    x, y, z = random.uniform(0,box_x), random.uniform(0,box_y), random.uniform(0,box_z)

    lines.append(f"{id} {id} 1 0 {x} {y} {z}")

out_file = os.path.join(script_dir, "water.data")
with open(out_file, "w") as f:
    f.write("\n".join(lines))

print("erfolgreich erstellt")



