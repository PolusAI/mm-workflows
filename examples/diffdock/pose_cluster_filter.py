# pylint: disable=E1101,E0401
"""Cluster predicted diffdock poses by centroid, \
    keep most confident from each cluster"""

import argparse
from pathlib import Path
import re

from rdkit import Chem
from rdkit.Chem import rdMolTransforms as rdmt
from rdkit.ML.Cluster import Butina
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('--predicted_poses', type=str, nargs='+', help='List of predicted pose files')
parser.add_argument('--centroid_cutoff', type=float,
                    help='cutoff for clusters of poses in same protein pocket')
args = parser.parse_args()

if args.predicted_poses is None:
    raise ValueError('No predicted poses provided')

predicted_poses = args.predicted_poses
# sanitize flag is still left True (default value), set removeHs to False to keep hydrogens
pred_mols = [Chem.SDMolSupplier(pose, removeHs=False)[0] for pose in predicted_poses]
if None in pred_mols:  # unreproducible failure rdkit returned None in virtual screen
    raise ValueError('Rdkit failed to read one of the poses')
predicted_poses = [Path(pose).name for pose in predicted_poses]


def EuclideanDist(pi: Chem.SDMolSupplier, pj: Chem.SDMolSupplier) -> float:
    """compute the Euclidean distance between two input molecule centroids

    Args:
        pi (Chem.SDMolSupplier): the ith molecule
        pj (Chem.SDMolSupplier): the jth molecule

    Returns:
        float: distance between two input molecule centroids
    """
    # GetConformers will just use input coordinates if conformations are not pre-generated
    confi = pi.GetConformers()[0]
    confj = pj.GetConformers()[0]
    centeri = rdmt.ComputeCentroid(confi)
    centerj = rdmt.ComputeCentroid(confj)
    dv = np.array([centeri.x, centeri.y, centeri.z]) - np.array([centerj.x, centerj.y, centerj.z])
    return float(np.sqrt(np.dot(dv, dv)))


def parse_confidence(file_name: str) -> float:
    """This function returns the confidence score from a filename.
    Filenames must follow the format 'rankX_confidenceY.mol',
    where X is a positive integer and Y is a float.
    This format is the default for DiffDock outputs."

    Args:
        file_name (str): The filename of output pose

    Returns:
        float: The confidence value from pose
    """
    confidence = float(re.findall('rank[0-9]+_confidence(.*).sdf', file_name)[0])
    return confidence


# Cluster a pose into group of other poses via centroid distance if beneath threshold
true_poses = []
index_to_cluster_index = {}
clusters = Butina.ClusterData(pred_mols, len(pred_mols), args.centroid_cutoff, distFunc=EuclideanDist)
for cluster_index, cluster in enumerate(clusters):
    names = [predicted_poses[index] for index in cluster]
    confidences = [parse_confidence(name) for name in names]
    max_confidence = max(confidences)
    confidence_to_indices = dict(zip(confidences, cluster))
    max_index = confidence_to_indices[max_confidence]
    true_poses.append(max_index)
    index_to_cluster_index[max_index] = cluster_index

# save pose index and cluster index to filtered_filename
lines = [f'{index} {index_to_cluster_index[index]} \n' for index in true_poses]
with open("filtered_poses.txt", 'w', encoding="utf-8") as file:
    file.writelines(lines)
