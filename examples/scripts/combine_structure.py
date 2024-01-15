import sys
import os
import argparse
from typing import Optional

from rdkit import Chem  # pylint: disable=import-error


def parse_arguments() -> argparse.Namespace:
    """ This function parses the arguments.

    Returns:
        argparse.Namespace: The command line arguments
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_structure1', required=True)
    parser.add_argument('--input_structure2', required=True)
    parser.add_argument('--output_structure_path', required=True)
    args = parser.parse_args()
    return args


def read_xyz_rdkit(input_structure_path: str) -> Optional[Chem.rdchem.Mol]:  # pylint: disable=c-extension-no-member
    """ Read a PDB file using RDKit

    Args:
        input_structure_path (str): The path to the xyz structure

    Returns:
        Optional[Chem.rdchem.Mol]: The created molecule object
    """
    xyz = Chem.rdmolfiles.MolFromXYZFile(input_structure_path)  # pylint: disable=c-extension-no-member

    if not xyz:
        print(f'Error: failed to generate molecule from file {input_structure_path}')
        return None

    return xyz


def combine_structure_rdkit(input_structure1_path: str, input_structure2_path: str, output_structure_path: str) -> None:
    """ Combine two structures into a single PDB file using RDKit

    Args:
        input_structure1_path (str): The path to the xyz structure 1
        input_structure2_path (str): The path to the xyz structure 2
        output_structure_path (str): The path to the output combined structure
    """

    structure1 = read_xyz_rdkit(input_structure1_path)
    structure2 = read_xyz_rdkit(input_structure2_path)

    if structure1 and structure2:
        combo = Chem.CombineMols(structure1, structure2)  # pylint: disable=no-member

        with Chem.PDBWriter(output_structure_path) as writer:  # pylint: disable=no-member
            writer.write(combo)


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
