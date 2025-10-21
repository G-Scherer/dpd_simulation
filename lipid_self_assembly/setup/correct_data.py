# This file corrects the molecule IDs in the "Atoms" section and the angle ID in the "Angles" section

import os

# Path to file
script_dir = os.path.dirname(os.path.abspath(__file__))

data_file = os.path.join(script_dir, "lammps_raw.data")

print(f"Bearbeite Datei: {data_file}")

# Read data
with open(data_file, 'r') as f:
    content = f.read()

# Split into list
lines = content.split('\n')
# Iniit new lines
new_lines = []

# Init flags
in_atoms_section = False
atoms_started = False
in_bonds_section = False
bonds_started = False
in_angles_section = False
angles_started = False

for line in lines:
    stripped = line.strip()
    #print(stripped)
    if stripped == "1  bond types":
        parts = line.split()
        parts[0] = str(2)
        print("bei bonds oben angekommen")

        new_lines.append("\t\t\t\t\t "+parts[0]+"  "+parts[1]+" "+parts[2])
        continue

    if stripped == "1  angle types":
        parts = stripped.split()
        parts[0] = str(3)

        new_lines.append("\t\t\t\t\t "+parts[0]+"  "+parts[1]+" "+parts[2])
        continue


    
    # Check for "Atoms" header
    if stripped == 'Atoms':
        in_atoms_section = True
        atoms_started = False
        new_lines.append(line)
        continue
    
    if in_atoms_section:
        # Skip empty line
        if not atoms_started and stripped == '':
            atoms_started = True
            new_lines.append(line)
            continue
        
        # Check for end of "Atoms" section
        if atoms_started and stripped == 'Bonds':
            in_atoms_section = False
            
        
        # Go into "Atoms" section
        if atoms_started and stripped:
            parts = line.split()
            if len(parts) >= 7:  
                atom_id = int(parts[0])
                
                # Calculate mol ID
                if atom_id <= 1200:
                    # POPC: 100 molecules, 12 atoms each
                    mol_id = ((atom_id - 1) // 12) + 1
                else:
                    # Water: begin after 100 POPC
                    mol_id = 100 + (atom_id - 1200)
                
                # Change mol ID 
                parts[1] = str(mol_id)
                
                # Save new lines
                new_lines.append(' '.join(parts))
                continue
    
    # Check for "Bonds" header
    if stripped == "Bonds":
        in_bonds_section = True
        bonds_started = False
        new_lines.append(line)
        continue

    if in_bonds_section:
        # Skip empty line
        if not bonds_started and stripped == "":
            bonds_started = True
            new_lines.append(line)
            continue

        if bonds_started and stripped == "Angles":
            in_bonds_section = False

        if bonds_started and stripped:
            parts = line.split()
            if len(parts) >= 4:
                bond_id = int(parts[0])
                # Starting at 3 every 11th bond type is 2

                if ((bond_id - 3) % 11 == 0):
                    parts[1] = str(2)

                new_lines.append(" ".join(parts))
                continue

    
    # Check for "Angles" header
    if stripped == "Angles":
        in_angles_section = True
        angles_started = False
        new_lines.append(line)
        continue

    if in_angles_section:
        # Skip empty line
        if not angles_started and stripped == "":
            angles_started = True
            new_lines.append(line)
            continue

        # Go into "Angles" section
        if angles_started and stripped:
            parts = line.split()
            if len(parts) >= 5:
                angle_id = int(parts[0])
                # Change angle type
                parts[1] = str(1)

                if ((angle_id - 1) % 8 == 0):
                    parts[1] = str(2)
                
                if ((angle_id - 5) % 8 == 0):
                    parts[1] = str(3)
                

                new_lines.append(" ".join(parts))
                continue


    # Leave other lines unchanged
    new_lines.append(line)

out_file = os.path.join(script_dir, "lammps.data")
# Write file
with open(out_file, 'w') as f:
    f.write('\n'.join(new_lines))

print(f"{data_file} wurde erfolgreich bearbeitet")
