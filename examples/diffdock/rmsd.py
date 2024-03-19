import argparse
from pathlib import Path
import json
from typing import Dict

from rdkit import Chem
from rdkit.Chem import AllChem

parser = argparse.ArgumentParser()
parser.add_argument('--predicted_poses', help='predicted pose')
parser.add_argument('--reference_pose', type=str, help='reference pose')
parser.add_argument('--output_json_name', type=str, help='output json file')
args = parser.parse_args()

predicted_poses = args.predicted_poses.split(',')

ref_mol = Chem.MolFromMolFile(args.reference_pose)
pred_mols = [Chem.MolFromMolFile(pose) for pose in predicted_poses]
predicted_poses = [Path(pose).name for pose in predicted_poses]


output_dict: Dict[str, Dict[str, float]] = {}
for i, pred_name in enumerate(predicted_poses):
    pred_mol = pred_mols[i]
    RMSD = AllChem.GetBestRMS(pred_mol, ref_mol)
    output_dict.update({pred_name: RMSD})

with open(args.output_json_name, 'w', encoding='utf-8') as f:
    json.dump(output_dict, f)
