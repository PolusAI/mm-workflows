# pylint: disable=no-member
import sys
import os
import argparse
from typing import List

import openmm.app as omma  # pylint: disable=import-error


def parse_arguments() -> argparse.Namespace:
    """ This function parses the arguments.

    Returns:
        argparse.Namespace: The command line arguments
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_pdb_path', type=str)
    parser.add_argument('--output_pdb_path', type=str)
    args = parser.parse_args()
    return args


def extract_protein(input_pdb_path: str, output_pdb_path: str) -> None:
    """ Extract protein from the PDB file

    Args:
        input_pdb_path (str): The path to the input pdb file
        output_pdb_path (str): The path to the output pdb file
    """

    pdb = omma.PDBFile(input_pdb_path)
    bonded_atom_idxs: List[int] = []
    for bond in pdb.topology.bonds():
        bonded_atom_idxs.extend([bond.atom1.index, bond.atom2.index])

    atom_idxs = [atom.index for atom in pdb.topology.atoms()]
    stray_atom_idxs = list(set(atom_idxs) - set(bonded_atom_idxs))

    stray_atoms: List[omma.topology.Atom] = []
    for atom in pdb.topology.atoms():
        if atom.index in stray_atom_idxs:
            stray_atoms.append(atom)

    modeller = omma.Modeller(pdb.topology, pdb.positions)
    modeller.delete(stray_atoms)

    # Delete water atoms
    modeller.deleteWater()

    pdb.topology = modeller.topology
    pdb.positions = modeller.positions

    with open(output_pdb_path, mode="w", encoding='utf-8') as wfile:
        omma.PDBFile.writeFile(pdb.topology, pdb.positions, wfile, keepIds=True)


def main() -> None:
    """ Reads the command line arguments and extract protein from the PDB file
    """
    args = parse_arguments()

    if not os.path.exists(args.input_pdb_path):
        print(f'Error: Can not find file {args.input_pdb_path}')
        sys.exit(1)

    extract_protein(args.input_pdb_path, args.output_pdb_path)


if __name__ == '__main__':
    main()
