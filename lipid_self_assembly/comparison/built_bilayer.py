import MDAnalysis as mda
import numpy as np
import freud
import matplotlib.pyplot as plt
import os

dir = "/home/gabriel/Dokumente/bachelor_thesis/dpd_simulation/lipid_self_assembly/comparison"

top_file = os.path.join(dir, "equilibration.gro")
traj_file = os.path.join(dir, "equilibration.xtc")

window_size = 3

u = mda.Universe(top_file, traj_file)

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

S_values = np.convolve(S_values, np.ones(window_size)/window_size, mode="valid")

plt.plot(S_values, color="black", linewidth=2, alpha=0.7)
plt.ylim((0,1))
plt.xlabel("frames")
plt.ylabel("Nematic order parameter")
plt.title("Lipid Orientation")

out_file = os.path.join(dir, "nematic_comp.png")

plt.savefig(out_file)
plt.show()