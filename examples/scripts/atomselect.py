from workflow_types import string, pdbfile  # pylint: disable=import-error
# NOTE: No other top-level imports supported


def main(selection_string: str, input_pdb_path: str, output_pdb_path: str) -> None:
    """Restrict a PDB file to a selection of atoms and save it.

    Args:
        selection_string (str): Selection string for mdtraj
        input_pdb_path (str): The path to the input PDB file
        output_pdb_path (str): The path to the output PDB file
    """
    import mdtraj  # pylint: disable=import-outside-toplevel
    traj = mdtraj.load(input_pdb_path)
    print(traj)
    selection_indices = traj.topology.select(selection_string)
    print(selection_indices)
    traj.restrict_atoms(selection_indices)
    traj.save(output_pdb_path)


inputs = {'selection_string': string,
          'input_pdb_path': pdbfile,
          'output_pdb_path': {**string, 'default': 'selection.pdb'}}
outputs = {'output_pdb_path': ('$(inputs.output_pdb_path)', pdbfile)}
