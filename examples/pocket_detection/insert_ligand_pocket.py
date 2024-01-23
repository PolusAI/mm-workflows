"""Insert ligand into protein pocket(s)"""

import argparse

from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np
from typing import List


parser = argparse.ArgumentParser()
parser.add_argument('--protein_filename', type=str, help='input protein file')
parser.add_argument('--ligand_filename', type=str, help='input ligand file')
parser.add_argument('--pocket_predictions', type=str, help='input pocket prediction csv file')
parser.add_argument('--top_n_pockets', type=str, help='top n pockets to insert ligand into')
args = parser.parse_args()

steric_threshold = .5  # Threshold for steric clashes
translation_refinement_increment = 0.5  # Increment for refining translation direction
max_translation_iter = 20  # Maximum number of iterations for refining translation direction
# the closer dot_prod_acceptable is to -1 the more the vector must point
# opposite of ligand center to protein center
dot_prod_acceptable = -0.05  # Acceptable dot product for solvent direction


def compute_bounding_box_volume(molecule: Chem.SDMolSupplier, confId: int, atom_ids: List[int]) -> float:
    """Compute the volume of the bounding box of a molecule

    Args:
        protein_lig (Chem.SDMolSupplier): The protein-ligand complex
        confId (int): Conformer ID
        atom_ids (List[int]): List of atom indices

    Returns:
        float: Volume of the bounding box
    """
    coords = np.array([molecule.GetConformer(confId).GetAtomPosition(i) for i in atom_ids])
    min_coords = np.min(coords, axis=0)
    max_coords = np.max(coords, axis=0)
    volume = float(np.prod(max_coords - min_coords))
    return volume


def check_steric_clash(protein_lig: Chem.SDMolSupplier, pro_indices: List[int], lig_indices: List[int], thresh: float) -> bool:
    """Check for steric clashes between protein and ligand atoms

    Args:
        protein_lig (Chem.SDMolSupplier): The protein-ligand complex
        pro_indices (List[int]): protein atom indices
        lig_indices (List[int]): ligand atom indices
        steric_threshold (float): Threshold for steric clashes

    Returns:
        bool: True if there is a steric clash, False otherwise
    """
    # compute pairwise distances between protein and ligand atoms
    distances = AllChem.Get3DDistanceMatrix(protein_lig, confId=0)

    protein_ligand_distances = distances[np.ix_(pro_indices, lig_indices)]

    # Check for steric clashes
    clash = bool(np.any(protein_ligand_distances < thresh))

    return clash


# assuming rank, center_x, center_y, center_z, surf_atom_ids are in csv file
rank_to_center = {}
rank_to_surf_atom_ids = {}
with open(args.pocket_predictions, 'r', encoding="utf-8") as f:
    lines = f.readlines()
    headers = lines[0].strip().split(',')
    headers = [x.strip() for x in headers]
    for line in lines[1:]:
        data = {headers[i]: value.strip() for i, value in enumerate(line.split(','))}
        rank = int(data['rank'])
        center_x = float(data['center_x'])
        center_y = float(data['center_y'])
        center_z = float(data['center_z'])
        surf_atom_ids = [int(x) - 1 for x in data['surf_atom_ids'].strip().split(' ')]
        rank_to_center[rank] = (center_x, center_y, center_z)
        rank_to_surf_atom_ids[rank] = surf_atom_ids

top_n = int(args.top_n_pockets)
top_n_centers = [rank_to_center[rank] for rank in range(1, top_n + 1)]

# rdkit valence error (sanitization error) is ignored
protein = Chem.MolFromPDBFile(args.protein_filename, removeHs=False, sanitize=False)
ligand = Chem.MolFromMolFile(args.ligand_filename, removeHs=False)

# Insert ligand into protein pocket(s)
# Compute bounding box of ligand and protein and compare volumes
# if ligand volume is larger than protein volume, then ligand cannot be inserted into protein
# Would be very complicated to do rigid rotation and conformer scanning,
# so we will just skip this case if can't insert via translation of center of ligand to center of pocket as first guess
# If protein pocket volume is larger but there is steric clash need to find translation direction
# that points towards solvent and translate ligand in that direction via tiny steps
# until sterically feasible (if possible)
for rank_idx, center in enumerate(top_n_centers):
    rank = rank_idx + 1
    ligand_copy = Chem.RWMol(ligand)
    conf = ligand_copy.GetConformer()

    # compute bounding box of ligand
    ligand_volume = compute_bounding_box_volume(ligand_copy, confId=0, atom_ids=list(range(ligand_copy.GetNumAtoms())))

    # define protein pocket coords as vectors from protein center to surf atom centers
    surf_atom_ids = rank_to_surf_atom_ids[rank]
    surf_atom_coords = np.array([protein.GetConformer(0).GetAtomPosition(i) for i in surf_atom_ids])
    surf_atom_coords = surf_atom_coords - np.array(center)

    # now compute bounding box of protein pocket
    protein_volume = compute_bounding_box_volume(protein, confId=0, atom_ids=surf_atom_ids)

    # compare volumes
    # note that alpha carbons are used as the surf atoms, so when translating towards boundary
    # can still clash with residue sticking into pocket
    print(f'Ligand volume {ligand_volume} protein volume {protein_volume}')
    if ligand_volume > protein_volume:
        continue
    else:
        print("Ligand should be able to fit into pocket")

    lig_num_atoms = conf.GetNumAtoms()
    # Now grab coordinates of ligand and determine ligand center position
    ligand_coords = np.array([conf.GetAtomPosition(i) for i in range(lig_num_atoms)])
    # Determine center of ligand
    ligand_center = [sum(x) / len(x) for x in zip(*ligand_coords)]
    # Determine translation vector
    translation_vector = [center[i] - ligand_center[i] for i in range(3)]
    # Translate ligand by using rdMolTransforms.TransformConformer
    new_coords = ligand_coords + translation_vector
    # Set the new coordinates back to the molecule
    for i in range(lig_num_atoms):
        conf.SetAtomPosition(i, new_coords[i])

    # Add ligand to protein
    pro_lig = Chem.CombineMols(protein, ligand_copy)
    # Get the number of atoms in the original protein and ligand
    num_protein_atoms = protein.GetNumAtoms()
    num_ligand_atoms = ligand_copy.GetNumAtoms()
    pro_idxs = list(range(num_protein_atoms))
    lig_idxs = list(range(num_protein_atoms, num_protein_atoms + num_ligand_atoms))
    steric_clash = check_steric_clash(pro_lig, pro_idxs, lig_idxs, steric_threshold)
    if steric_clash:
        print(f'Rank {rank} has steric clash, attempting to pull towards solvent')
        # Grab list of vectors from pocket center that point to all pocket protein atoms
        vectors = np.array([protein.GetConformer(0).GetAtomPosition(i) - np.array(center) for i in surf_atom_ids])
        # Convert to unit vectors
        unit_vectors = vectors / np.linalg.norm(vectors, axis=0)

        # find center of entire protein
        pro_center = np.array([sum(x) / len(x)
                              for x in zip(*[protein.GetConformer(0).GetAtomPosition(i) for i in pro_idxs])])

        lig_center_to_pro_center = pro_center - np.array(center)
        # normalize the vector
        lig_center_to_pro_center = lig_center_to_pro_center / np.linalg.norm(lig_center_to_pro_center)
        # compute dot products of all unit_vectors to lig_center_to_pro_center
        dot_products = np.dot(unit_vectors, lig_center_to_pro_center)
        # if the dot product is negative, then the vector points away from protein center.
        # this assumes that the pocket is not near the center of the protein
        solvent_vectors = unit_vectors[dot_products < dot_prod_acceptable]
        # for troubleshooting print the proteins atoms that point towards solvent
        # solvent_protein_atoms = [surf_atom_ids[i] for i in range(len(surf_atom_ids)) \
        # if dot_products[i] < dot_prod_acceptable]

        # find the vector that points the most towards solvent
        # average the vectors that point towards solvent
        solvent_vector = np.mean(solvent_vectors, axis=0)
        # normalize the vector
        solvent_vector = solvent_vector / np.linalg.norm(solvent_vector)

        # Iteratively translate ligand in translation direction via tiny steps until sterically feasible
        iterations = 0
        conf = pro_lig.GetConformer()
        new_coords = np.array([conf.GetAtomPosition(i) for i in lig_idxs])
        Chem.MolToPDBFile(pro_lig, f'initial_trans{rank}.pdb', confId=0)
        while steric_clash:
            if iterations > max_translation_iter:
                break
            new_coords = new_coords + solvent_vector * translation_refinement_increment
            # only update ligand coordinates (ligand is after protein indices)
            for i in range(num_protein_atoms, num_protein_atoms + num_ligand_atoms):
                lig_index = i - num_protein_atoms
                conf.SetAtomPosition(i, new_coords[lig_index])
            steric_clash = check_steric_clash(pro_lig, pro_idxs, lig_idxs, steric_threshold)
            iterations += 1
        if not steric_clash:
            Chem.MolToPDBFile(pro_lig, f'protein_with_ligand_rank_{rank}.pdb', confId=0)

    else:
        Chem.MolToPDBFile(pro_lig, f'protein_with_ligand_rank_{rank}.pdb', confId=0)
