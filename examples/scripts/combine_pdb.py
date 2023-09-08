# pylint: disable=no-member
import sys
import os
import argparse
from typing import Optional
import numpy as np

from rdkit import Chem
import mdtraj as mdj
import parmed as pmd

def parse_arguments() -> argparse.Namespace:
    """ This function parses the arguments.

    Returns:
        argparse.Namespace: The command line arguments
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_structure1', required=True)
    parser.add_argument('--input_structure2', required=True)
    parser.add_argument('--output_structure_path', required=True)
    parser.add_argument('--method', required=False, default='mdtraj')
    args = parser.parse_args()
    return args

def read_pdb_rdkit(input_structure_path: str) -> Optional[Chem.rdchem.Mol]:
    """ Read a PDB file using RDKit 

    Args:
        input_structure_path (str): The path to the pdb structure 

    Returns:
        Optional[Chem.rdchem.Mol]: The created molecule object
    """

    pdb = Chem.MolFromPDBFile(input_structure_path, sanitize=False, removeHs=False, proximityBonding=False)

    if not pdb:
        print(f'Error: failed to generate molecule from file {input_structure_path}')
        sys.exit(1)

    return pdb

def combine_pdb_rdkit(input_structure1_path: str, input_structure2_path: str, output_structure_path: str) -> None:
    """ Combine two PDB structures into a single PDB file using RDKit

    Args:
        input_structure1_path (str): The path to the pdb structure 1 
        input_structure2_path (str): The path to the pdb structure 2 
        output_structure_path (str): The path to the output combined structure
    """

    pdb1 = read_pdb_rdkit(input_structure1_path)
    pdb2 = read_pdb_rdkit(input_structure2_path)

    combo = Chem.CombineMols(pdb1, pdb2)

    with Chem.PDBWriter(output_structure_path) as writer:
        writer.write(combo)

def combine_pdb_mdtraj(input_structure1_path: str, input_structure2_path: str, output_structure_path: str) -> None:
    """ Combine two PDB structures into a single PDB file using MDtraj

    Args:
        input_structure1_path (str): The path to the pdb structure 1 
        input_structure2_path (str): The path to the pdb structure 2 
        output_structure_path (str): The path to the output combined structure
    """

    pdb1 = mdj.load_pdb(input_structure1_path)
    pdb2 = mdj.load_pdb(input_structure2_path)

    top1, xyz1 = pdb1.top, pdb1.xyz
    top2, xyz2 = pdb2.top, pdb2.xyz

    top1_2 = top1.join(top2)
    xyz1_2 = np.concatenate((xyz1, xyz2), axis=1)

    combined_pdb = mdj.Trajectory(xyz=xyz1_2, topology=top1_2)
    combined_pdb.save_pdb(output_structure_path)

def combine_pdb_parmed(input_structure1_path: str, input_structure2_path: str, output_structure_path: str) -> None:
    """ Combine two PDB structures into a single PDB file using Parmed

    Args:
        input_structure1_path (str): The path to the pdb structure 1 
        input_structure2_path (str): The path to the pdb structure 2 
        output_structure_path (str): The path to the output combined structure
    """

    # The output of the load_file function, as explained in here
    # https://parmed.github.io/ParmEd/html/readwrite.html
    # is contingent on the input file type. To effectively combine two structures,
    # the load_file function should return a "Structure" object.
    # For instance, when working with Mol2/Mol3 files, the structure argument should
    # be set to True to obtain a Structure object.
    pdb1 = pmd.load_file(input_structure1_path, structure=True)
    pdb2 = pmd.load_file(input_structure2_path, structure=True)

    # The overloaded "+" operator here only functions properly when
    # the input objects are instances of the "Structure" object
    # https://parmed.github.io/ParmEd/html/structure.html#structure-combining
    combined = pdb1 + pdb2

    combined.save(output_structure_path, format='pdb', overwrite=True)

def main() -> None:
    """ Reads the command line arguments and combine two PDB structures into a single PDB file
    """
    args = parse_arguments()

    if not os.path.exists(args.input_structure1):
        print(f'Error: Can not find file {args.input_structure1}')
        sys.exit(1)

    if not os.path.exists(args.input_structure2):
        print(f'Error: Can not find file {args.input_structure2}')
        sys.exit(1)

    if args.method == 'mdtraj':
        combine_pdb_mdtraj(args.input_structure1, args.input_structure2, args.output_structure_path)
    elif args.method == 'rdkit':
        combine_pdb_rdkit(args.input_structure1, args.input_structure2, args.output_structure_path)
    elif args.method == 'parmed':
        combine_pdb_parmed(args.input_structure1, args.input_structure2, args.output_structure_path)
    else:
        print('The method must be one of the following: mdtraj, rdkit, or parmed.')

if __name__ == '__main__':
    main()
