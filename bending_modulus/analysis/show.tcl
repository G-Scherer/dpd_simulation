# Load topology and trajectory
mol new membrane.gro type gro waitfor all
mol addfile production.xtc type xtc waitfor all

mol delrep 0 top

# White Background
color Display Background white
display projection Orthographic
axes location Off
pbc box

# Define custom colors
color change rgb 0 0.2 0.6 0.8    ;# blue for head groups
color change rgb 1 0.8 0.3 0.2    ;# orange for tails

# Head groups representation
mol representation VDW 0.5 30
mol color ColorID 0
mol material AOShiny
mol selection "resname POPC and name NC3 PO4 GL1 GL2"
mol addrep top

# Tail representation
mol representation VDW 0.3 30
mol color ColorID 1
mol material AOShiny
mol selection "resname POPC and name C1A C2A C3A C4A D2A C5A C1B C2B C3B C4B D2B C5B"
mol addrep top