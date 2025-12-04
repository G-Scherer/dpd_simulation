# Load trajectories
mol new trajectories_256.lammpstrj type lammpstrj waitfor all

mol delrep 0 top

# White Background
color Display Background white
display projection Perspective
axes location Off


# Color and Shape
mol representation VDW 0.1 30
mol color Type
mol material AOShiny
mol addrep top

color Type 4 23
color Type 1 10
color Type 2 6
color Type 3 2
color Type 5 30
color Type 6 3
