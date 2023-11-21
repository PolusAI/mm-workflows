"""Ranks predicted poses from DiffDock via confidence score"""
import argparse
import re
from typing import Dict

parser = argparse.ArgumentParser()
parser.add_argument('--diffdock_poses', type=str, nargs='+', help='List of poses from DiffDock')
parser.add_argument('--top_n_confident', type=int, help='Top n confident poses')
parser.add_argument('--top_percent_confidence', type=float, help='top confidence percent cutoff')
args = parser.parse_args()


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


confidences = [parse_confidence(name) for name in args.diffdock_poses]
confidence_to_pose: Dict[float, str] = dict(zip(confidences, args.diffdock_poses))
file_to_index: Dict[str, int] = dict(zip(args.diffdock_poses, range(len(args.diffdock_poses))))

# First filter by absolute value top_n_confident
# if user only wants to use top_percent_confidence, can set top_n_confident to trivially high number
# if user only wants to use top_n_confident, then can set top_percent_confidence to 100
sorted_list = sorted(confidence_to_pose.items(), reverse=True)
sorted_pose_list = [v for (k, v) in sorted_list]
poses = sorted_pose_list[:args.top_n_confident]
# Next filter by top percentage of confident poses
num_poses = int(args.top_percent_confidence*.01*len(poses))
poses = poses[:num_poses]
# Write to ranked_filename
with open("ranked_poses.txt", 'w', encoding="utf-8") as file:
    file.writelines([f"{file_to_index[string]} \n" for string in poses])
