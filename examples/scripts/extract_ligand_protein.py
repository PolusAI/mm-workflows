# pylint: disable=no-member
import sys
import os
import argparse

import MDAnalysis as mda


def parse_arguments() -> argparse.Namespace:
    """ This function parses the arguments.

    Returns:
        argparse.Namespace: The command line arguments
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_pdb_path', type=str)
    parser.add_argument('--output_pdb_path', type=str)
    parser.add_argument('--output_pdb_ligand_path', type=str)
    args = parser.parse_args()
    return args


def extract_ligand_protein(input_pdb_path: str, output_pdb_path: str, output_pdb_ligand_path: str) -> None:
    """ Extract ligand & protein from the PDB file

    Args:
        input_pdb_path (str): The path to the input pdb file
        output_pdb_path (str): The path to the output pdb file
        output_pdb_ligand_path (str): The path to the output pdb ligand file
    """

    # Load the PDB file
    u = mda.Universe(input_pdb_path)

    # Get unique residue names
    protein_atoms = u.select_atoms('protein')  # use simple atom selection when possible

    # Create a new Universe with only protein atoms
    protein_u = mda.Universe.empty(n_atoms=protein_atoms.n_atoms, trajectory=True)  # needed for coordinates
    protein_u.atoms = protein_atoms

    # duplicate the universe object
    dup_u = mda.Universe(input_pdb_path)

    # now do the same for the ligand, not protein and not water or salts
    ligand_atoms = u.select_atoms('not protein')

    try:
        # guess the bonds, since input PDB may not have bonds
        dup_u.atoms.guess_bonds()
    except ValueError:
        # ValueError: vdw radii for types: AS. These can be defined manually using the keyword 'vdwradii'
        print('Error: Could not guess bonds. Check the input PDB file.')

    has_bonds = False
    try:
        num_bonds = len(dup_u.atoms.bonds)
        has_bonds = True
    except mda.exceptions.NoDataError:
        print('No bonds found in the PDB file.')

    # Identify water molecules based on the connectivity pattern (Oxygen bonded to two Hydrogens)
    if has_bonds:
        water_indices = set()
        for atom in dup_u.atoms:  # dont use selection resname == 'HOH', pdb file may have different water residue names
            if atom.name == 'O' and len(atom.bonds) == 2:  # if hydrogens are added
                bonded_atoms_names = set([a.name for a in atom.bonded_atoms])
                if bonded_atoms_names == {'H'}:  # Check if both bonds are Hydrogens
                    water_indices.add(atom.index)
                    water_indices.update([a.index for a in atom.bonded_atoms])

        # now want to remove all salts, waters without H
        non_bonded = set()
        for atom in dup_u.atoms:
            if len(atom.bonds) == 0:
                non_bonded.add(atom.index)

        # Remove water by excluding the water indices
        if len(water_indices) > 0:
            water_indices_string = ' '.join([str(i) for i in water_indices])
            ligand_atoms = ligand_atoms.select_atoms(f'not index {water_indices_string}')

        # Remove non bonded atoms
        if len(non_bonded) > 0:
            non_bonded_string = ' '.join([str(i) for i in non_bonded])
            ligand_atoms = ligand_atoms.select_atoms(f'not index {non_bonded_string}')

    ligand_u = mda.Universe.empty(n_atoms=ligand_atoms.n_atoms, trajectory=True)  # needed for coordinates
    ligand_u.atoms = ligand_atoms

    with open(output_pdb_path, mode="w", encoding='utf-8') as wfile:
        protein_u.atoms.write(output_pdb_path)
    if len(ligand_u.atoms) > 0:  # will crash if no ligand atoms
        with open(output_pdb_ligand_path, mode="w", encoding='utf-8') as wfile:
            ligand_u.atoms.write(output_pdb_ligand_path)


def main() -> None:
    """ Reads the command line arguments and extract protein from the PDB file
    """
    args = parse_arguments()

    if not os.path.exists(args.input_pdb_path):
        print(f'Error: Can not find file {args.input_pdb_path}')
        sys.exit(1)

    extract_ligand_protein(args.input_pdb_path, args.output_pdb_path, args.output_pdb_ligand_path)


if __name__ == '__main__':
    main()
