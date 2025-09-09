import sys
import os
import argparse
import ast
from typing import Dict


def parse_arguments() -> argparse.Namespace:
    """ This function parses the arguments.

    Returns:
        argparse.Namespace: The command line arguments
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('--input_path', type=str)
    parser.add_argument('--output_path', type=str)
    parser.add_argument('--residue_name', type=str)
    parser.add_argument('--column_index', type=int, required=False, default=None)
    parser.add_argument('--line_start_column_idxs', type=str, required=False, default=None)
    args = parser.parse_args()
    return args


def rename_residues(input_path: str, output_path: str, residue_name: str,
                    column_index: int, line_start_column_idxs: str) -> None:
    """ Modifies the names of residues based on the specified column index or 
    line start and column indices

    Args:
        input_path (str): input file path
        output_path (str): output file path
        residue_name (str): The new residue name to which all residues should be changed 
        column_index (int): The index of column
        line_start_column_idxs (str): A dictionary string indicating the line start and corresponding column indices, 
        such as '{"ATOM": 3, "HETATM": 3}'.
    """
    line_start_column_idxs_dict: Dict[str, int] = {}

    if line_start_column_idxs:
        # Parse the dictionary string to a dictionary
        line_start_column_idxs_dict = ast.literal_eval(line_start_column_idxs)

    with open(input_path, mode='r', encoding='utf-8') as f:
        lines = f.readlines()

    lines_new = []
    for line in lines:
        l = line
        words = l.split()
        if column_index:
            if len(words) >= column_index:
                l = l.replace(words[column_index], residue_name)
        elif line_start_column_idxs:
            for line_start, cindex in line_start_column_idxs_dict.items():
                if words[0] == line_start:
                    l = l.replace(words[cindex], residue_name)

        lines_new.append(l)

    with open(output_path, mode='w', encoding='utf-8') as f:
        f.writelines(lines_new)


def main() -> None:
    """ Reads the command line arguments and renames the residues
    """
    args = parse_arguments()

    if args.input_path and not os.path.exists(args.input_path):
        print(f'Error: Can not find file {args.input_path}')
        sys.exit(1)

    if (args.column_index is None) and (args.line_start_column_idxs is None):
        print("Error: No column index information is provided")
        sys.exit(1)

    if (args.column_index) and (args.line_start_column_idxs):
        print("Error: One of column_index or line_start_column_idxs arguments should be provided")
        sys.exit(1)

    rename_residues(args.input_path, args.output_path, args.residue_name,
                    args.column_index, args.line_start_column_idxs)


if __name__ == '__main__':
    main()
