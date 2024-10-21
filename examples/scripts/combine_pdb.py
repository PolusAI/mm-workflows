import sys
import os
import argparse
from typing import Optional
import numpy as np

import MDAnalysis as mda


def parse_arguments() -> argparse.Namespace:
    """ This function parses the arguments.

    Returns:
        argparse.Namespace: The command line arguments
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_structure1', type=str)
    parser.add_argument('--input_structure2', type=str)
    parser.add_argument('--output_structure_path', type=str)
    args = parser.parse_args()
    return args


def combine_structure_rdkit(input_structure1_path: str, input_structure2_path: str, output_structure_path: str) -> None:
    """ Combine two structures into a single PDB file using MDAnalysis

    Args:
        input_structure1_path (str): The path to the PDB structure 1 
        input_structure2_path (str): The path to the PDB structure 2 
        output_structure_path (str): The path to the output combined PDB structure
    """

    u1 = mda.Universe(input_structure1_path)
    u2 = mda.Universe(input_structure2_path)

    merged = mda.Merge(u1.atoms, u2.atoms)
    merged_elements = np.concatenate([u1.atoms.elements, u2.atoms.elements])

    merged.atoms.elements = merged_elements
    all_atoms = merged.select_atoms("all")
    all_atoms.write(output_structure_path)


def main() -> None:
    """ Reads the command line arguments and combine two structures into a single PDB file
    """
    args = parse_arguments()

    if not os.path.exists(args.input_structure1):
        print(f'Error: Can not find file {args.input_structure1}')
        sys.exit(1)

    if not os.path.exists(args.input_structure2):
        print(f'Error: Can not find file {args.input_structure2}')
        sys.exit(1)

    combine_structure_rdkit(args.input_structure1, args.input_structure2, args.output_structure_path)


if __name__ == '__main__':
    main()
