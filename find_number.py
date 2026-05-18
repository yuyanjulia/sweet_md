def find_residue_change():
    with open("brazzein_processed.gro", 'r') as file:
        lines = file.readlines()[1:]

    receptor_start = 1
    receptor_end = None
    ligand_start = None
    ligand_end = None

    for line in lines:
        parts = line.split()
        if len(parts) == 1:
            ligand_end = int(parts[0])

        if len(parts) == 5:
            residue = int(parts[0][:-3])

            if residue == 1:
                atom_number = int(parts[1][-5:])
                ligand_start = atom_number
                receptor_end = atom_number - 1
                break

    with open("atom_number.txt", 'w') as output_file:
        output_file.write(f"{receptor_start}\n")
        output_file.write(f"{receptor_end}\n")
        output_file.write(f"{ligand_start}\n")
        output_file.write(f"{ligand_end}\n")

find_residue_change()
