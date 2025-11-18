# Load trajectories
mol new trajectories.lammpstrj type lammpstrj waitfor all

mol delrep 0 top

# White Background
color Display Background white
display projection Perspective
axes location Off


# Color and Shape
mol representation VDW 0.1 12
mol color Type
mol material Opaque
mol addrep top

color Type 4 23
color Type 1 16
color Type 2 28
color Type 3 6
color Type 5 1
color Type 6 18
