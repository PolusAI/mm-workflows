"""Generates pymol session files with pdb and ligand"""
import argparse
from pathlib import Path

from pymol import cmd

parser = argparse.ArgumentParser()
parser.add_argument('--ligand_path', type=str, help='Ligand file path')
parser.add_argument('--pdb_path', type=str, help='protein PDB file path')
args = parser.parse_args()

pdb_name = Path(args.pdb_path).name
pdb_name = pdb_name.split('.')[0]
cmd.load(args.pdb_path, 'pdb')

if args.ligand_path is not None:
    ligand_name = Path(args.ligand_path).name
    ligand_name = ligand_name.split('.')[0]
    cmd.load(args.ligand_path)
    cmd.show('sticks', ligand_name)
    cmd.show('sticks', 'all within 5 of '+ligand_name)
    cmd.label('n. CA within 5 of '+ligand_name, 'resn')
    cmd.zoom(ligand_name)
    name = f'{pdb_name}_{ligand_name}'
else:
    name = pdb_name

cmd.save(f'{name}.pse')
