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
all_frames_box = []

for frame in u.trajectory:

    frame_mol_positions = []

    for mol in u.residues:

        # only lipids
        if len(mol.atoms) == 12:
            positions_nm = mol.atoms.positions.copy() / 10.0
            frame_mol_positions.append(positions_nm)

    all_frames_mol_positions.append(frame_mol_positions)
    all_frames_box.append(u.dimensions[:3].copy())

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


def com(mol, box):

    box = np.asarray(box, dtype=float)

    # choose arbitrary reference bead
    ref = mol[0].copy()
    
    # calculate distance
    delta = mol - ref

    # apply nearest image convention, if delta/box > 0.5 then we subtract a box length for that bead
    # everything is then in intervall [-box/2, box/2]
    delta = delta - box * np.round(delta/box)
    unwrapped = ref + delta 

    return unwrapped.mean(axis=0)

# calulate contact

# ...existing code...

# calulate contact

contact_number = []

for frame in all_frames_mol_positions:

    # calculate com of every mol
    box = np.array(u.dimensions[:3])
    coms = np.array([com(mol, box) for mol in frame])

    # matrix of pairwise distance vecotrs
    delta = coms[:,None,:] - coms[None,:,:]

    # apply pbc
    delta = delta - box*np.round(delta/box)

    # take norm and discard distance with self
    dists = np.linalg.norm(delta, axis=2)
    np.fill_diagonal(dists, np.inf)

    # discard dists bigger than cutoff and get number of close neighbors for each lipid
    neighbors_per_lipid = (dists <= 1.1).sum(axis=1)

    # get mean of every lipid
    contact_number.append(neighbors_per_lipid.mean())

contact_number = np.array(contact_number)

print(f"\n=== Contact Number Summary ===")
print(f"Total frames: {len(contact_number)}")
print(f"Contact number stats:")
print(f"  Mean: {contact_number.mean():.4f}")
print(f"  Std:  {contact_number.std():.4f}")
print(f"  Min:  {contact_number.min():.4f}")
print(f"  Max:  {contact_number.max():.4f}")
print(f"  First 10 values: {contact_number[:10]}")

# ...existing code (convolution etc.)...

S_values = np.convolve(S_values, np.ones(window_size)/window_size, mode="valid")
contact_number = np.convolve(contact_number, np.ones(window_size)/window_size, mode="valid")


plt.plot(S_values, color="black", linewidth=2, alpha=0.7)
plt.ylim((0,1))
plt.xlabel("frames")
plt.ylabel("Nematic order parameter")
plt.title("Lipid Orientation")

out_file = os.path.join(dir, "nematic_comp.png")

plt.savefig(out_file)
plt.show()

plt.plot(contact_number, color="black", linewidth=2, alpha=0.7)
plt.xlabel("frames")
plt.ylabel("Contact number")
plt.title("Contact number")

out_file = os.path.join(dir, "contact_number.png")

plt.savefig(out_file)
plt.show()