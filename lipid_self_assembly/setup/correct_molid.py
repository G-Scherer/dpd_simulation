import os

script_dir = os.path.dirname(os.path.abspath(__file__))
data_file = os.path.join(script_dir, "lammps.data")

print(f"Bearbeite Datei: {data_file}")

# Lese Datei
with open(data_file, 'r') as f:
    content = f.read()

# Teile in Zeilen
lines = content.split('\n')
# Neue Zeilen
new_lines = []
in_atoms_section = False
atoms_started = False

for line in lines:
    stripped = line.strip()
    #print(stripped)
    
    # Prüfe ob "Atoms" Header kommt
    if stripped == 'Atoms':
        print(line, "true gesetzt")
        in_atoms_section = True
        atoms_started = False
        new_lines.append(line)
        continue
    
    # Wenn wir in Atoms-Sektion sind
    if in_atoms_section:
        # Leerzeile nach "Atoms" überspringen
        if not atoms_started and stripped == '':
            atoms_started = True
            #new_lines.append(line)
            continue
        
        # Ende der Atoms-Sektion (Leerzeile oder neue Sektion)
        if atoms_started and stripped == 'Bonds':
            print(line, "false gesetzt")
            in_atoms_section = False
            new_lines.append(line)
            continue
        
        # Verarbeite Atom-Zeile
        if atoms_started and stripped:
            parts = line.split()
            if len(parts) >= 7:  # Gültige Atom-Zeile
                atom_id = int(parts[0])
                print("ist bei >7 angekommen")
                
                # Berechne Molekül-ID
                if atom_id <= 1536:
                    # POPC: 128 Moleküle, je 12 Atome
                    mol_id = ((atom_id - 1) // 12) + 1
                else:
                    # Wasser: nach den 128 POPC-Molekülen
                    mol_id = 128 + (atom_id - 1536)
                
                # Ersetze Molekül-ID (zweites Element)
                parts[1] = str(mol_id)
                
                # Schreibe neue Zeile
                new_lines.append(' '.join(parts))
                continue
    
    # Alle anderen Zeilen unverändert
    new_lines.append(line)

# Schreibe zurück
with open(data_file, 'w') as f:
    f.write('\n'.join(new_lines))
