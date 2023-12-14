import sys
import os
import argparse
from typing import List, Tuple

from pdbfixer import PDBFixer
import openmm.app as omma


def parse_arguments() -> argparse.Namespace:
    """ This function parses the arguments.

    Returns:
        argparse.Namespace: The command line arguments
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_pdb_path', required=False, type=str, default=None)
    parser.add_argument('--output_pdb_path', type=str)
    parser.add_argument('--input_helper_pdb_path', required=False, type=str, default=None)
    parser.add_argument('--add_atoms', type=str, choices=('all', 'heavy', 'none'))
    parser.add_argument('--add_residues', required=False, action='store_true')
    parser.add_argument('--pdbid', required=False, type=str, default=None)
    parser.add_argument('--url', required=False, type=str, default=None)
    parser.add_argument('--replace_nonstandard', required=False, action='store_true')
    parser.add_argument('--keep_heterogens', required=False, type=str, default='all', choices=('all', 'water', 'none'))
    args = parser.parse_args()
    return args


def find_missing_residues(fixer: PDBFixer) -> PDBFixer:
    """ Finds the missing residues and adds missing residues within a 
    chain to prevent "floppy tails," which can lead to an increase in the box size, 
    significantly increasingthe computation time. This step is taken as floppy tails 
    are generally not critical for binding.

    Args:
        fixer (PDBFixer): The input PDBFixer instance

    Returns:
        PDBFixer: The output PDBFixer instance with added missing residues
    """
    fixer.findMissingResidues()
    fixer_chains: List[omma.topology.Chain] = list(fixer.topology.chains())
    fixer_chain_res_idx_pairs: List[Tuple[int, int]] = list(fixer.missingResidues.keys())
    for chain_idx, res_idx in fixer_chain_res_idx_pairs:
        chain = fixer_chains[chain_idx]
        if res_idx == 0 or res_idx == len(list(chain.residues())):
            del fixer.missingResidues[tuple([chain_idx, res_idx])]
    return fixer


def runpdbfixer(input_pdb_path: str, input_helper_pdb_path: str, output_pdb_path: str,
                add_atoms: str, add_res: bool, pdbid: str, url: str, rep_nonstandard: bool, heterogens: str) -> None:
    """ Fixes the protein structure using PDBFixer.PDBFixer offers options 
    to add hydrogens and solvate the system, but in our usage, we employ 
    PDBFixer solely for adding missing heavy atoms and residues.

    Args:
        input_pdb_path (str): The input PDB structure path
        output_pdb_path (str): The output PDB structure path
        input_helper_pdb_path (str): The input helper PDB structure path
        add_atoms (str): What missing atoms to add: all, heavy, hydrogen, or none
        add_res (bool): If set to True, adds missing residues
        pdbid (str): PDB id from RCSB 
        url (str): URL to retrieve PDB from
        rep_nonstandard (bool): Replace nonstandard residues with standard equivalents
    """
    # The input can be one of these options
    if input_pdb_path:
        fixer = PDBFixer(filename=input_pdb_path)
    elif pdbid:
        fixer = PDBFixer(pdbid=pdbid)
    elif url:
        fixer = PDBFixer(url=url)

    if add_res:
        if input_helper_pdb_path:
            helper_fixer = PDBFixer(filename=input_helper_pdb_path)
            # Finds the missing residues based on the PDBbind structure but uses
            # the sequence info from the helper file
            helper_fixer.topology = fixer.topology
            helper_fixer = find_missing_residues(helper_fixer)
            fixer.missingResidues = helper_fixer.missingResidues
        else:
            fixer = find_missing_residues(fixer)
    else:
        fixer.missingResidues = {}

    if rep_nonstandard:
        fixer.findNonstandardResidues()
        fixer.replaceNonstandardResidues()

    if heterogens == 'none':
        fixer.removeHeterogens(False)
    elif heterogens == 'water':
        fixer.removeHeterogens(True)

    fixer.findMissingAtoms()
    if add_atoms not in ('all', 'heavy'):
        fixer.missingAtoms = {}
        fixer.missingTerminals = {}

    # Adds identified missing atoms and residues
    fixer.addMissingAtoms()
    with open(output_pdb_path,  mode="w", encoding='utf-8') as wfile:
        omma.PDBFile.writeFile(fixer.topology, fixer.positions, wfile, keepIds=True)


def main() -> None:
    """ Reads the command line arguments and combine two structures into a single PDB file
    """
    args = parse_arguments()

    if args.input_pdb_path and not os.path.exists(args.input_pdb_path):
        print(f'Error: Can not find file {args.input_pdb_path}')
        sys.exit(1)

    if args.input_helper_pdb_path and not os.path.exists(args.input_helper_pdb_path):
        print(f'Error: Can not find file {args.input_helper_pdb_path}')
        sys.exit(1)

    if (args.input_pdb_path is None) and (args.pdbid is None) and (args.url is None):
        print("Error: No input is provided")
        sys.exit(1)
    runpdbfixer(args.input_pdb_path, args.input_helper_pdb_path,
                args.output_pdb_path, args.add_atoms, args.add_residues,
                args.pdbid, args.url, args.replace_nonstandard, args.keep_heterogens)


if __name__ == '__main__':
    main()
