import argparse
import os
import numpy as np
import random
from scipy.constants import Boltzmann,Avogadro

parser = argparse.ArgumentParser(description='Create lipid soup system')

# define arguments, box default taken from past sim
parser.add_argument('--lipids', type=int, default=5000, help='Number of lipids')
parser.add_argument('--density', type=float, default=3.0, help='Number density of beads')
parser.add_argument("--apl", type=float, default=0.879, help="apl in nm^2")
parser.add_argument("--inner_radius", type=float, default=10, help="inner radius in nm")

# parse arguments
args = parser.parse_args()

# save variables
n_lipids = args.lipids
density = args.density
apl = args.apl
R_input = args.inner_radius

r_ref = 0.711 #nm
e_ref = Boltzmann*298.15 #K
q = 8.861242189860825 

apl = apl / r_ref**2

R_in = R_input/r_ref
R_out = R_in + 2.8/r_ref
R = (R_in + R_out)/2

box_z = (n_lipids*apl)/(4*np.pi*R)

N_out = int((2*np.pi*R_out*box_z)/apl)
N_in = n_lipids - N_out

box_x = box_y = 2*R_out + 12

print(f"radius: R={R}, innen Radius R_in={R_in}, außen Radius {R_out} \n lipids {n_lipids}, außen {N_out}, innnen {N_in}, box {box_y, box_x,box_z}")

box_vol = box_x*box_y*box_z

n_water = int(box_vol*density - 12*n_lipids)

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

n_z = int(box_z/np.sqrt(apl))
n_theta_out = int(np.ceil(N_out/n_z))
n_theta_in = int(np.ceil(N_in/n_z))

print(f"n_z {n_z}, n_theta_out {n_theta_out}, n_theta_in {n_theta_in}")

def get_radial_rotation(theta, inward=False):
    # Lipid von Z in XY-Ebene kippen (90 Grad um Y)
    R_tilt = np.array([[0,0,1], [0,1,0], [-1,0,0]])
    # Ausrichtung am Radiusvektor
    phi = theta + np.pi if inward else theta
    R_z = np.array([[np.cos(phi), -np.sin(phi), 0],
                    [np.sin(phi),  np.cos(phi), 0],
                    [0, 0, 1]])
    return R_z.dot(R_tilt)

mol_id = 1
center_x, center_y = box_x / 2, box_y / 2

for leaflet in ['in', 'out']:
    current_R = R_in if leaflet == 'in' else R_out
    current_n_theta = n_theta_in if leaflet == 'in' else n_theta_out
    is_inward = (leaflet == 'in')
    
    for i in range(n_z):
        for j in range(current_n_theta):
            if (leaflet == 'in' and mol_id > N_in) or (leaflet == 'out' and mol_id > n_lipids):
                break
                
            z_com = (i + 0.5) * (box_z / n_z) + random.uniform(-0.1, 0.1)
            # Staggering (hexagonalere Packung)
            shift = (np.pi / current_n_theta) if i % 2 == 0 else 0
            theta = j * (2 * np.pi / current_n_theta) + shift + random.uniform(-0.02, 0.02)
            
            x_com = center_x + current_R * np.cos(theta)
            y_com = center_y + current_R * np.sin(theta)
            b = np.array([x_com, y_com, z_com])
            
            R_mat = get_radial_rotation(theta, inward=is_inward)
            
            for n, coords in enumerate(popc_centered):
                c = R_mat.dot(coords) + b
                # Z-PBC sicherstellen
                c[2] = c[2] % box_z
                lines.append(f"{(mol_id-1)*12+n+1} {mol_id} {bead_types[n]} {charges[n]} {c[0]} {c[1]} {c[2]}")
            mol_id += 1

water_added = 0
atom_id = n_lipids * 12 + 1
water_mol_id = n_lipids + 1
margin = 1.2 / r_ref # Sicherheitsabstand zur Membran

while water_added < n_water:
    x, y, z = random.uniform(0, box_x), random.uniform(0, box_y), random.uniform(0, box_z)
    dist = np.sqrt((x - center_x)**2 + (y - center_y)**2)
    
    # Wasser nur innen (r < R_in) oder außen (r > R_out)
    if dist < (R_in - margin) or dist > (R_out + margin):
        lines.append(f"{atom_id} {water_mol_id} {bead_types[-1]} {charges[-1]} {x} {y} {z}")
        water_added += 1
        atom_id += 1
        water_mol_id += 1

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

out_file = os.path.join(script_dir, f"cylinder_{R_input}.data")
with open(out_file, "w") as f:
    f.write("\n".join(lines))

print("erfolgreich erstellt")



