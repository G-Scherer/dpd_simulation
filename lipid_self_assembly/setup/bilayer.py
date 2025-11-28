import argparse
import os
import numpy as np
import random
from scipy.constants import Boltzmann,Avogadro

parser = argparse.ArgumentParser(description='Create lipid soup system')

# define arguments, box default taken from past sim
parser.add_argument('--lipids', type=int, default=256, help='Number of lipids')
parser.add_argument('--density', type=float, default=3.0, help='Number density of beads')
parser.add_argument("--length", type=float, default=10.6, help="box length in nm")

# parse arguments
args = parser.parse_args()

# save variables
n_lipids = args.lipids
density = args.density
box = args.length

r_ref = 0.711 #nm
e_ref = Boltzmann*298.15 #K
q = 8.861242189860825 

box_x = box_y = box/r_ref

box_z = box/r_ref*1.5

box_vol = box_x*box_y*box_z

n_water = int(box_vol*density - 12*n_lipids)


# ...existing code...

random.seed(75)

print(f"Creating system:")
print(f"  Lipids: {n_lipids}")
print(f"  Water: {n_water}")
print(f"  Box volume: {box_vol} reduced units")
print(f" density: {(n_lipids*12 + n_water)/box_vol}")

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

lines.append(f"0 {box_x} xlo xhi")
lines.append(f"0 {box_y} ylo yhi")
lines.append(f"0 {box_z} zlo zhi")

lines.append("")

lines.append("Masses")

lines.append("")

for i in range(1,7):
    lines.append(f"{i} 1")

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
popc_mol /= r_ref
popc_centered = popc_mol - np.mean(popc_mol, axis = 0)

# bead types[n] for n=0,...,11 in popc and n=-1 is water
bead_types = [5, 6, 3, 3, 1, 2, 1, 1, 1, 1, 1, 1, 4]
charges = np.array([1.0, -1.0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
charges *= q

# ✅ BILAYER SETUP
# Berechne wie viele Lipide pro Leaflet (grid)
n_per_leaflet = n_lipids // 2
grid_size = int(np.ceil(np.sqrt(n_per_leaflet)))

# Grid spacing in x,y
spacing_x = box_x / grid_size
spacing_y = box_y / grid_size

z_center = box_z * 0.5  # Mitte der Box
leaflet_separation = 3.0  # Abstand zwischen Leaflets in reduced units (~2.1 nm)

z_lower = z_center - leaflet_separation / 2  # Unteres Leaflet
z_upper = z_center + leaflet_separation / 2  # Oberes L

print(f"Creating bilayer structure:")
print(f"  Grid: {grid_size} x {grid_size}")
print(f"  Leaflet spacing: {spacing_x:.2f} x {spacing_y:.2f}")
print(f"  Lower leaflet z: {z_lower:.2f}")
print(f"  Upper leaflet z: {z_upper:.2f}")

mol_id = 1

# ✅ Unteres Leaflet (Heads zeigen nach unten, -z Richtung)
for i in range(grid_size):
    for j in range(grid_size):
        if mol_id > n_per_leaflet:
            break
        
        # Position im Grid mit kleiner Störung
        x_com = (i + 0.5) * spacing_x + random.uniform(-0.3, 0.3)
        y_com = (j + 0.5) * spacing_y + random.uniform(-0.3, 0.3)
        z_com = z_lower + random.uniform(-0.2, 0.2)
        
        # Periodische Randbedingungen
        x_com = x_com % box_x
        y_com = y_com % box_y
        
        b = np.array([x_com, y_com, z_com])
        
        # Rotation: Heads zeigen nach unten (-z)
        # Lipid ist ~entlang z-Achse, drehe um 180° in x
        theta_x = np.pi + random.uniform(-0.2, 0.2)  # ~180° ± noise
        theta_y = random.uniform(-0.2, 0.2)
        theta_z = random.uniform(0, 2*np.pi)
        
        R_x = np.array([[1,0,0],[0,np.cos(theta_x),-np.sin(theta_x)],[0,np.sin(theta_x),np.cos(theta_x)]])
        R_y = np.array([[np.cos(theta_y),0,-np.sin(theta_y)],[0,1,0],[np.sin(theta_y),0,np.cos(theta_y)]])
        R_z = np.array([[np.cos(theta_z),-np.sin(theta_z),0],[np.sin(theta_z),np.cos(theta_z),0],[0,0,1]])
        R = R_x.dot(R_y).dot(R_z)
        
        for n, coords in enumerate(popc_centered):
            coords = R.dot(coords)
            coords = coords + b
            coords[0] = coords[0] % box_x
            coords[1] = coords[1] % box_y
            coords[2] = coords[2] % box_z
            
            lines.append(f"{(mol_id - 1)*12+n+1} {mol_id} {bead_types[n]} {charges[n]} {coords[0]} {coords[1]} {coords[2]}")
        
        mol_id += 1

# ✅ Oberes Leaflet (Heads zeigen nach oben, +z Richtung)
for i in range(grid_size):
    for j in range(grid_size):
        if mol_id > n_lipids:
            break
        
        # Position im Grid mit kleiner Störung (leicht versetzt zum unteren Leaflet)
        x_com = (i + 0.5) * spacing_x + random.uniform(-0.3, 0.3)
        y_com = (j + 0.5) * spacing_y + random.uniform(-0.3, 0.3)
        z_com = z_upper + random.uniform(-0.2, 0.2)
        
        x_com = x_com % box_x
        y_com = y_com % box_y
        
        b = np.array([x_com, y_com, z_com])
        
        # Rotation: Heads zeigen nach oben (+z)
        # Lipid bleibt ~entlang z-Achse
        theta_x = random.uniform(-0.2, 0.2)  # ~0° ± noise
        theta_y = random.uniform(-0.2, 0.2)
        theta_z = random.uniform(0, 2*np.pi)
        
        R_x = np.array([[1,0,0],[0,np.cos(theta_x),-np.sin(theta_x)],[0,np.sin(theta_x),np.cos(theta_x)]])
        R_y = np.array([[np.cos(theta_y),0,-np.sin(theta_y)],[0,1,0],[np.sin(theta_y),0,np.cos(theta_y)]])
        R_z = np.array([[np.cos(theta_z),-np.sin(theta_z),0],[np.sin(theta_z),np.cos(theta_z),0],[0,0,1]])
        R = R_x.dot(R_y).dot(R_z)
        
        for n, coords in enumerate(popc_centered):
            coords = R.dot(coords)
            coords = coords + b
            coords[0] = coords[0] % box_x
            coords[1] = coords[1] % box_y
            coords[2] = coords[2] % box_z
            
            lines.append(f"{(mol_id - 1)*12+n+1} {mol_id} {bead_types[n]} {charges[n]} {coords[0]} {coords[1]} {coords[2]}")
        
        mol_id += 1

# ✅ WASSER: Fülle Bereiche außerhalb der Bilayer
water_mol_id_start = n_lipids + 1
atom_id = n_lipids * 12 + 1

# Wasser oben (z > z_upper + margin)
water_added = 0
margin = 2.0  # Abstand zur Bilayer

while water_added < n_water // 2 and water_added < n_water:
    x_com = random.uniform(0, box_x)
    y_com = random.uniform(0, box_y)
    z_com = random.uniform(z_upper + margin, box_z)
    
    coords = [x_com, y_com, z_com]
    
    lines.append(f"{atom_id} {water_mol_id_start + water_added} {bead_types[-1]} {charges[-1]} {coords[0]} {coords[1]} {coords[2]}")
    
    atom_id += 1
    water_added += 1

# Wasser unten (z < z_lower - margin)
while water_added < n_water:
    x_com = random.uniform(0, box_x)
    y_com = random.uniform(0, box_y)
    z_com = random.uniform(0, z_lower - margin)
    
    coords = [x_com, y_com, z_com]
    
    lines.append(f"{atom_id} {water_mol_id_start + water_added} {bead_types[-1]} {charges[-1]} {coords[0]} {coords[1]} {coords[2]}")
    
    atom_id += 1
    water_added += 1

lines.append("")

# ...existing code (Bonds, Angles bleiben gleich)...

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

lines.append("Bond Coeffs")

lines.append("")

k_b = 1250 * 1/2 * 1000 * 1/Avogadro * 1/e_ref * r_ref**2

lines.append(f"1 {k_b} {0.47/r_ref}")
lines.append(f"2 {k_b} {0.37/r_ref}")

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

lines.append("")

lines.append("Angle Coeffs")

lines.append("")

k_a1 = 25 * 1/2 * 1000 * 1/Avogadro * 1/e_ref 
k_a2 = 45 * 1/2 * 1000 * 1/Avogadro * 1/e_ref 

lines.append(f"1 {k_a1} 180")
lines.append(f"2 {k_a1} 120")
lines.append(f"3 {k_a2} 120")

out_file = os.path.join(script_dir, "bilayer.data")
with open(out_file, "w") as f:
    f.write("\n".join(lines))

print("erfolgreich erstellt")



