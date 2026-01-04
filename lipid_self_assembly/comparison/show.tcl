# Load trajectories
mol new comp_128_high_T.lammpstrj type lammpstrj waitfor all

mol delrep 0 top

# White Background
color Display Background white
display projection Orthographic
axes location Off
pbc box


# Color and Shape
mol representation VDW 0.1 30
mol color Type
mol material AOShiny
mol selection "not type 4 and not type 5 and not type 6"
mol addrep top

mol representation VDW 0.2 30
mol color Type
mol material AOShiny
mol selection "type 5 or type 6"
mol addrep top

color Name 4 23
color Name 1 10
color Name 2 6
color Name 3 2
color Name 5 30
color Name 6 3
