import os

def get_element(atom_name, res_name):
    """
    Returns the correct element based on atom name and residue name.
    Based on standard PDB nomenclature for amino acids.
    """
    # Strip any leading/trailing spaces from atom name
    atom_name = atom_name.strip()
    
    # Standard backbone atoms
    backbone = {
        'N': 'N',
        'CA': 'C',
        'C': 'C',
        'O': 'O',
        'OXT': 'O'
    }
    
    # Common sidechain atoms
    sidechain = {
        'CB': 'C',
        'CG': 'C',
        'CG1': 'C',
        'CG2': 'C',
        'CD': 'C',
        'CD1': 'C',
        'CD2': 'C',
        'CE': 'C',
        'CE1': 'C',
        'CE2': 'C',
        'CE3': 'C',
        'CZ': 'C',
        'CZ2': 'C',
        'CZ3': 'C',
        'CH2': 'C',
        'ND1': 'N',
        'ND2': 'N',
        'NE': 'N',
        'NE1': 'N',
        'NE2': 'N',
        'NH1': 'N',
        'NH2': 'N',
        'NZ': 'N',
        'OD1': 'O',
        'OD2': 'O',
        'OE1': 'O',
        'OE2': 'O',
        'OG': 'O',
        'OG1': 'O',
        'OH': 'O',
        'SD': 'S',
        'SG': 'S',
        'SE': 'Se'
    }
    
    # Special cases for specific residues
    special_cases = {
        'HIS': {
            'HD1': 'H',
            'HE2': 'H'
        },
        'CYS': {
            'SG': 'S'
        },
        'MET': {
            'SD': 'S'
        },
        'MSE': {  # Selenomethionine
            'SE': 'Se'
        }
    }
    
    # Check backbone atoms first
    if atom_name in backbone:
        return backbone[atom_name]
    
    # Check residue-specific special cases
    if res_name in special_cases and atom_name in special_cases[res_name]:
        return special_cases[res_name][atom_name]
    
    # Check common sidechain atoms
    if atom_name in sidechain:
        return sidechain[atom_name]
    
    # Handle hydrogen atoms
    if atom_name.startswith('H'):
        return 'H'
    
    # For any other carbon atoms starting with C
    if atom_name.startswith('C'):
        return 'C'
    
    # For any other nitrogen atoms starting with N
    if atom_name.startswith('N'):
        return 'N'
    
    # For any other oxygen atoms starting with O
    if atom_name.startswith('O'):
        return 'O'
    
    # Default case
    return atom_name[0]

def fix_pdb_file(input_file, output_file):
    """Fixes a single PDB file by adding missing element entries."""
    fixed_lines = []
    modified = False
    
    with open(input_file, 'r') as infile:
        for line in infile:
            if line.startswith("ATOM") or line.startswith("HETATM"):
                atom_name = line[12:16].strip()
                res_name = line[17:20].strip()
                element = get_element(atom_name, res_name)
                
                # Check if the element column is already filled
                if line[76:78].strip() != element:
                    modified = True
                    # Ensure the line is at least 76 characters long by padding with spaces
                    if len(line) < 76:
                        line = line.rstrip() + " " * (76 - len(line.rstrip()))
                    # Add element in columns 77-78
                    new_line = line[:76] + f"{element:>2}"
                    # Preserve any remaining content after column 78
                    if len(line) > 78:
                        new_line += line[78:]
                    if not new_line.endswith('\n'):
                        new_line += '\n'
                    fixed_lines.append(new_line.rstrip('\n'))
                else:
                    fixed_lines.append(line.rstrip('\n'))
            else:
                fixed_lines.append(line.rstrip('\n'))
    
    # Write to output file if modified
    if modified:
        with open(output_file, 'w') as outfile:
            outfile.write('\n'.join(fixed_lines) + '\n')
    return modified

def process_pdb_folder(input_folder, output_folder):
    """Processes all PDB files in a folder, fixing missing element entries."""
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    
    pdb_files = [f for f in os.listdir(input_folder) if f.endswith('.pdb')]
    print(f"Found {len(pdb_files)} PDB file(s) in folder '{input_folder}'.")
    
    for pdb_file in pdb_files:
        input_path = os.path.join(input_folder, pdb_file)
        output_path = os.path.join(output_folder, pdb_file)
        
        print(f"Processing '{pdb_file}'...", end=" ")
        modified = fix_pdb_file(input_path, output_path)
        
        if modified:
            print("fixed and saved.")
        else:
            print("no changes needed.")
    print(f"Fixed PDB files are saved to '{output_folder}'.")

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="Fix missing element entries in PDB files.")
    parser.add_argument("input_folder", help="Path to the folder containing PDB files.")
    parser.add_argument("output_folder", help="Path to the folder to save fixed PDB files.")
    
    args = parser.parse_args()
    process_pdb_folder(args.input_folder, args.output_folder)