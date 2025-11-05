import MDAnalysis as mda
import os
import numpy as np
import freud
import matplotlib.pyplot as plt
import argparse

parser = argparse.ArgumentParser(description="Calculate nematic order")

parser.add_argument("--equil", type=str, required=True, help="filepath for equilibriation from lipid_self_assembly")
parser.add_argument("--prod", type=str, required=True, help="filepath for production from lipid_self_assembly")
parser.add_argument("--top", type=str, required=True, help="filepath from lipid_self_assembly")
parser.add_argument("--cutoff", type=float, required=True, help="cutoff value for contact number")

args = parser.parse_args()

top = args.top
equil = args.equil
prod = args.prod
cutoff = args.cutoff

# window size for rolling average
window_size = 10


dir = "/home/gabriel/Dokumente/bachelor_thesis/dpd_simulation/lipid_self_assembly/"
top_file = os.path.join(dir, top)
equil_file = os.path.join(dir, equil)
prod_file = os.path.join(dir, prod)

u = mda.Universe(top_file, format="DATA", atom_style = "id resid type charge x y z")

print("=== Topology loaded ===")
print(f"Total atoms: {len(u.atoms)}")
print(f"Total residues: {len(u.residues)}")

u.load_new(equil_file, format="LAMMPSDUMP", dt=1.0)

equil_frames = len(u.trajectory)

all_frames_mol_positions = []

for frame in u.trajectory:

    frame_mol_positions = []

    for mol in u.residues:

        # only lipids
        if len(mol.atoms) == 12:
            frame_mol_positions.append(mol.atoms.positions.copy())

    all_frames_mol_positions.append(frame_mol_positions)

u.load_new(prod_file, format="LAMMPSDUMP", dt=1.0)

prod_frames = len(u.trajectory)

for frame in u.trajectory:

    frame_mol_positions = []

    for mol in u.residues:

        # only lipids
        if len(mol.atoms) == 12:
            frame_mol_positions.append(mol.atoms.positions.copy())

    all_frames_mol_positions.append(frame_mol_positions)

# [frame id] [molid] [atom id in mol] [axis] as selection for all_frames_mol_positions

print("info shape: ", np.shape(all_frames_mol_positions))

n_lipids = np.shape(all_frames_mol_positions)[1]

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
    neighbors_per_lipid = (dists <= cutoff).sum(axis=1)

    # get mean of every lipid
    contact_number.append(neighbors_per_lipid.mean())

contact_number = np.array(contact_number)

# calculate rolling average

S_values = np.convolve(S_values, np.ones(window_size)/window_size, mode="valid")
contact_number = np.convolve(contact_number, np.ones(window_size)/window_size, mode="valid")

middle = int(len(S_values)/2)

plt.plot(np.arange(middle),S_values[:middle], linewidth=2, alpha=0.7, label="equilibriation")
plt.plot(np.arange(middle-1,2*middle),S_values[middle:], linewidth=2, alpha=0.7, label="production")
plt.ylim((0,1))
plt.xlabel("frames")
plt.ylabel("Nematic order parameter")
plt.title("Lipid Orientation")

plt.grid()
plt.legend()
plt.savefig("nematic_order.png")
plt.show()

plt.plot(np.arange(middle),contact_number[:middle], linewidth=2, alpha=0.7, label="equilibriation")
plt.plot(np.arange(middle-1, 2*middle),contact_number[middle:], linewidth=2, alpha=0.7, label="production")
plt.xlabel("frames")
plt.ylabel("Contact number")
plt.title("Contact number")

plt.grid()
plt.legend()
plt.savefig("contact_number.png")
plt.show()


    
    


    





