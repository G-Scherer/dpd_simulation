# Load trajectories
mol new membrane.gro type gro waitfor all
mol addfile production.xtc type xtc waitfor all

mol delrep 0 top

# White Background
color Display Background white
display projection Orthographic
axes location Off
pbc box

color Name W 23
color Name C1A 10
color Name C3A 10
color Name C4A 10
color Name C1B 10
color Name C2B 10
color Name C3B 10
color Name C4B 10
color Name D2A 6
color Name GL1 2
color Name GL2 2
color Name NC3 30
color Name PO4 3


# Color and Shape
mol representation VDW 0.7 30
mol color ColorID 10
mol material AOShiny
mol selection "name C1A C3A C4A C1B C2B C3B C4B"
mol addrep top

mol representation VDW 0.7 30
mol color ColorID 6
mol material AOShiny
mol selection "name D2A"
mol addrep top

mol representation VDW 0.7 30
mol color ColorID 2
mol material AOShiny
mol selection "name GL1 GL2"
mol addrep top

mol representation VDW 1.4 30
mol color ColorID 3
mol material AOShiny
mol selection "name PO4"
mol addrep top

mol representation VDW 1.4 30
mol color ColorID 30
mol material AOShiny
mol selection "name NC3"
mol addrep top


