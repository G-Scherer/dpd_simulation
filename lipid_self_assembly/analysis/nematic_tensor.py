import MDAnalysis as mda
import os
import numpy as np
import freud
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser(description="Calculate nematic order")

parser.add_argument("--traj", type=str, required=True, help="filepath from lipid_self_assembly")
parser.add_argument("--top", type=str, required=True, help="filepath from lipid_self_assembly")

args = parser.parse_args()

top = args.top
traj = args.traj

dir = "/home/gabriel/Dokumente/bachelor_thesis/dpd_simulation/lipid_self_assembly/"
top_file = os.path.join(dir, top)
traj_file = os.path.join(dir, traj)

u = mda.Universe(top_file, format="DATA", atom_style = "id resid type charge x y z")

print("=== Topology loaded ===")
print(f"Total atoms: {len(u.atoms)}")
print(f"Total residues: {len(u.residues)}")

u.load_new(traj_file, format="LAMMPSDUMP", dt=1.0)

print("\n=== Trajectory loaded ===")
print(f"Total frames: {len(u.trajectory)}")
print(f"Box dimensions: {u.dimensions}")

all_frames_mol_positions = []

for frame in u.trajectory:

    frame_mol_positions = []

    for mol in u.residues:

        # only lipids
        if len(mol.atoms) == 12:
            frame_mol_positions.append(mol.atoms.positions.copy())

    all_frames_mol_positions.append(frame_mol_positions)

# [frame id] [molid] [atom id in mol] [axis] as selection for all_frames_mol_positions

# calculate allignement vector 

def u_vector(mol):

    head = mol[0]
    tail = (mol[7] + mol[11])/2

    u = (tail - head)/np.linalg.norm(tail - head)

    return u

nematic = freud.order.Nematic()

S_values = []
directors = []

for i, frame in enumerate(all_frames_mol_positions):

    orientations = np.array([u_vector(mol) for mol in frame])

    nematic.compute(orientations)

    S_values.append(nematic.order)
    directors.append(nematic.director)

    # if i % 10 == 0:
    #     print(f"frame {i}: S = {nematic.order}")

S_values = np.array(S_values)
directors = np.array(directors)

plt.plot(S_values, color="black", linewidth=2, alpha=0.7)
plt.ylim((0,1))
plt.xlabel("frames")
plt.ylabel("Nematic order parameter")
plt.title("Lipid Orientation")

traj_basename = os.path.splitext(os.path.basename(traj))[0]
output_filename = f"nematic_order_{traj_basename}.png"

plt.savefig(output_filename)
plt.show()


    
    


    





