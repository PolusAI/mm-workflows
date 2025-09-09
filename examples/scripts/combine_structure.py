import sys
import os
import zipfile
import argparse
from typing import Optional

from rdkit import Chem
import MDAnalysis as mda


GROMACS_INCLUDE_PATH = os.path.realpath('/opt/conda/share/gromacs/top')
TEMP_PATH = os.path.realpath('./')


def parse_arguments() -> argparse.Namespace:
    """ This function parses the arguments.

    Returns:
        argparse.Namespace: The command line arguments
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_structure1', required=True)
    parser.add_argument('--input_structure2', required=True)
    parser.add_argument('--input_top_zip_path', required=True)
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


def fix_atom_names(input_top_zip_path: str, output_structure_path: str) -> None:
    """ Modifies the atom names based on the topology file

    Args:
        input_top_zip_path (str): The path to the input zip file
        output_structure_path (str): The path to the output combined structure
    """

    with zipfile.ZipFile(input_top_zip_path, "r") as zip_read:
        zip_read.extractall(path=TEMP_PATH)

    file_list = [str(os.path.join(TEMP_PATH, f)) for f in zip_read.namelist()]
    top_file = next(name for name in file_list if name.endswith(".top"))

    u1 = mda.Universe(top_file, topology_format='ITP', include_dir=GROMACS_INCLUDE_PATH)
    u2 = mda.Universe(os.path.join(TEMP_PATH, 'temp_out.pdb'))

    u2.atoms.names = u1.atoms.names
    all_atoms = u2.select_atoms("all")
    all_atoms.write(output_structure_path)


def combine_structure_rdkit(input_structure1_path: str, input_structure2_path: str) -> None:
    """ Combine two structures into a single PDB file using RDKit and saves it in the temporary PDB file

    Args:
        input_structure1_path (str): The path to the xyz structure 1 
        input_structure2_path (str): The path to the xyz structure 2 
    """

    structure1 = read_xyz_rdkit(input_structure1_path)
    structure2 = read_xyz_rdkit(input_structure2_path)

    if structure1 and structure2:
        combo = Chem.CombineMols(structure1, structure2)  # pylint: disable=no-member

        with Chem.PDBWriter(os.path.join(TEMP_PATH, 'temp_out.pdb')) as writer:
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

    combine_structure_rdkit(args.input_structure1, args.input_structure2)
    fix_atom_names(args.input_top_zip_path, args.output_structure_path)


if __name__ == '__main__':
    main()
