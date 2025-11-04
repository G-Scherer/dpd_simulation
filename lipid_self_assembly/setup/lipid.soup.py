import argparse
import os
import numpy as np
import random

parser = argparse.ArgumentParser(description='Create lipid soup system')

# define arguments
parser.add_argument('--lipids', type=int, default=100, help='Number of lipids')
parser.add_argument('--water', type=int, default=1000, help='Number of water beads')
parser.add_argument('--density', type=float, default=3.0, help='Number density of beads')

# parse arguments
args = parser.parse_args()

# save variables
n_lipids = args.lipids
n_water = args.water
density = args.density
box_vol = (n_lipids*12+n_water)/density
box_length = box_vol**(1/3)

random.seed(123)

print(f"Creating system:")
print(f"  Lipids: {n_lipids}")
print(f"  Water: {n_water}")
print(f"  Box length: {box_length} reduced units")

script_dir = os.path.dirname(os.path.abspath(__file__))

lines = []

# Add header and basic information about system

lines.append("Lipid soup setup")

lines.append("")

lines.append(f"{12*n_lipids + n_water} atoms")
lines.append(f"{11*n_lipids} bonds")
lines.append(f"{8*n_lipids} angles")
lines.append("0 dihedrals")
lines.append("0 impropers")

lines.append("")

lines.append("6 atom types")
lines.append("2 bond types")
lines.append("3 angle types")
lines.append("0 dihedral types")
lines.append("0 improper types")

lines.append("")

for ax in ["x","y","z"]:
    lines.append(f"{0} {box_length} {ax}lo {ax}hi")

lines.append("")

lines.append("Masses")

lines.append("")

for i in range(1,7):
    lines.append(f"{i} 72.000000")

lines.append("")

lines.append("Atoms")

lines.append("")

# molecule setup from .gro, popc_mol[n,m] is nth(0,...,11) bead in mol, m is x,y,z(0,1,2)
popc_mol = np.array([[0.561, 0.534, 2.457],
                     [0.435, 0.553, 2.166],
                     [0.478, 0.512, 1.813],
                     [0.778, 0.470, 1.805],
                     [0.415, 0.469, 1.422],
                     [0.451, 0.703, 1.201],
                     [0.412, 0.577, 0.903],
                     [0.415, 0.606, 0.559],
                     [0.938, 0.611, 1.512],
                     [0.919, 0.457, 1.196],
                     [0.924, 0.591, 0.907],
                     [0.927, 0.557, 0.559]])

# divide by r_ref=0.711 nm to make unitless
popc_mol /= 0.711
popc_centered = popc_mol - np.mean(popc_mol, axis = 0)

# bead types[n] for n=0,...,11 in popc and n=-1 is water
bead_types = [5, 6, 3, 3, 1, 2, 1, 1, 1, 1, 1, 1, 4]
charges = [1.0, -1.0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

for mol_id in range(1,(n_lipids+1)):

    x_com, y_com, z_com = random.uniform(0,box_length), random.uniform(0,box_length), random.uniform(0,box_length)
    b = np.array([x_com, y_com, z_com])

    theta_x, theta_y, theta_z = random.uniform(0, 2*np.pi), random.uniform(0, 2*np.pi), random.uniform(0, 2*np.pi)
    R_x = np.array([[1,0,0],[0,np.cos(theta_x),-np.sin(theta_x)],[0,np.sin(theta_x),np.cos(theta_x)]])
    R_y = np.array([[np.cos(theta_y),0,-np.sin(theta_y)],[0,1,0],[np.sin(theta_y),0,np.cos(theta_y)]])
    R_z = np.array([[np.cos(theta_z),-np.sin(theta_z),0],[np.sin(theta_z),np.cos(theta_z),0],[0,0,1]])
    R = R_x.dot(R_y).dot(R_z)
    #mol_coords = popc_mol + np.array([x_com, y_com, z_com])
    #mol_coords = R.dot(mol_coords)
    #mol_coords = mol_coords % box_length

    for n,coords in enumerate(popc_centered):
        
        coords = R.dot(coords)
        coords = coords+b
        coords = coords % box_length

        lines.append(f"{(mol_id - 1)*12+n+1} {mol_id} {bead_types[n]} {charges[n]} {coords[0]} {coords[1]} {coords[2]}") 


for mol_id,atom_id in enumerate(range(n_lipids*12+1, 12*n_lipids+n_water+1),n_lipids):

    x_com, y_com, z_com = random.uniform(0,box_length), random.uniform(0,box_length), random.uniform(0,box_length)
    coords = [x_com, y_com, z_com]

    lines.append(f"{atom_id} {mol_id+1} {bead_types[-1]} {charges[-1]} {coords[0]} {coords[1]} {coords[2]}")

lines.append("")

lines.append("Bonds")

lines.append("")



for mol_id in range(n_lipids):

    for bond in range(1,12):

        if bond == 3:

            lines.append(f"{11*mol_id+bond} 2 {12*mol_id+bond} {12*mol_id+bond+1}")
            continue
        
        if bond == 4:

            lines.append(f"{11*mol_id+bond} 1 {12*mol_id+bond-1} {12*mol_id+bond+1}")
            continue
        
        if bond == 8:

            lines.append(f"{11*mol_id+bond} 1 {12*mol_id+bond-4} {12*mol_id+bond+1}")
            continue

        else:

            lines.append(f"{11*mol_id+bond} 1 {12*mol_id+bond} {12*mol_id+bond+1}")
        
lines.append("")

lines.append("Angles")

lines.append("")

for mol_id in range(n_lipids):

    for angle in range(1,9):

        if angle == 1:

            lines.append(f"{8*mol_id+angle} 2 {12*mol_id+angle+1} {12*mol_id+angle+2} {12*mol_id+angle+3}")
        
        if angle == 2:

            lines.append(f"{8*mol_id+angle} 1 {12*mol_id+angle} {12*mol_id+angle+1} {12*mol_id+angle+3}")

        if angle == 3:
            
            lines.append(f"{8*mol_id+angle} 1 {12*mol_id+angle} {12*mol_id+angle+2} {12*mol_id+angle+3}")
        
        if angle == 4:

            lines.append(f"{8*mol_id+angle} 3 {12*mol_id+angle+1} {12*mol_id+angle+2} {12*mol_id+angle+3}")
        
        if angle == 5:

            lines.append(f"{8*mol_id+angle} 1 {12*mol_id+angle+1} {12*mol_id+angle+2} {12*mol_id+angle+3}")

        if angle == 6:

            lines.append(f"{8*mol_id+angle} 1 {12*mol_id+angle-2} {12*mol_id+angle+3} {12*mol_id+angle+4}")

        if angle == 7:

            lines.append(f"{8*mol_id+angle} 1 {12*mol_id+angle+2} {12*mol_id+angle+3} {12*mol_id+angle+4}")

        if angle == 8:

            lines.append(f"{8*mol_id+angle} 1 {12*mol_id+angle+2} {12*mol_id+angle+3} {12*mol_id+angle+4}")



out_file = os.path.join(script_dir, "lammps.data")
with open(out_file, "w") as f:
    f.write("\n".join(lines))

print("erfolgreich erstellt")



