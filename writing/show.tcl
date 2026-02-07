# Load trajectories
mol new popc_viz.lammpstrj type lammpstrj waitfor all

mol delrep 0 top

# White Background
color Display Background white
display projection Orthographic
axes location Off
pbc box

# Color and Shape
# Atom ID 1-4: rot
mol representation VDW 0.1 30
mol color ColorID 1
mol material AOShiny
mol selection "serial 1 to 4"
mol addrep top

# Atom ID 5,6,9,10: lila
mol representation VDW 0.1 30
mol color ColorID 11
mol material AOShiny
mol selection "serial 5 6 9 10"
mol addrep top

# Atom ID 7,8,11,12: blau
mol representation VDW 0.1 30
mol color ColorID 0
mol material AOShiny
mol selection "serial 7 8 11 12"
mol addrep top
