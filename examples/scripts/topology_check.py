import argparse
import signal
from rdkit import Chem
from rdkit.Chem.rdchem import Atom, Mol


def map_atoms(mol1: Mol, mol2: Mol) -> dict:
    """
    Map atoms from the first molecule to corresponding atoms in the second molecule.

    Args:
        mol1 (Mol): RDKit molecule from the first PDB file.
        mol2 (Mol): RDKit molecule from the second PDB file.

    Returns:
        dict: Mapping of atoms from mol1 to mol2.
    """
    mol1_num_atoms = mol1.GetNumAtoms()
    mol2_num_atoms = mol2.GetNumAtoms()
    # if the mol being mapped has less atom then a match won't be found
    if mol1_num_atoms > mol2_num_atoms:
        match_atoms = mol1.GetSubstructMatch(mol2)
        mol2_index_to_mol1_index = {mol2_index: mol1_index for mol2_index, mol1_index in enumerate(match_atoms)}
        mol1_index_to_mol2_index = {k: v for v, k in mol2_index_to_mol1_index.items()}
    else:
        match_atoms = mol2.GetSubstructMatch(mol1)
        mol1_index_to_mol2_index = {mol1_index: mol2_index for mol1_index, mol2_index in enumerate(match_atoms)}
    if len(match_atoms) == 0:
        raise ValueError("No atoms matched between the two molecules.")

    return mol1_index_to_mol2_index


def compare_atom_properties(atom1: Atom, atom2: Atom) -> dict:
    """
    Compare properties of two atoms and return a dictionary of differences.

    Args:
        atom1 (Atom): RDKit atom from the first molecule.
        atom2 (Atom): RDKit atom from the second molecule.

    Returns:
        dict: Dictionary of property differences between the two atoms.
    """
    properties_diff = {}
    atom1_properties = get_properties(atom1)
    atom2_properties = get_properties(atom2)
    for key, atom1_value in atom1_properties.items():
        if atom1_value != atom2_properties[key]:
            properties_diff[key] = (atom1_value, atom2_properties[key])
    return properties_diff


def compare_atoms(mol1: Mol, mol2: Mol, atom_map: dict, intended_changes: list = []) -> dict:
    """
    Compare properties of mapped atoms in two molecules and return a dictionary of differences.

    Args:
        mol1 (Mol): RDKit molecule from the first PDB file.
        mol2 (Mol): RDKit molecule from the second PDB file.
        atom_map (dict): Mapping of atoms from mol1 to mol2.
        intended_changes (list): List of intended changes (e.g., "added_hydrogen").

    Returns:
        dict: Dictionary of property differences between mapped atoms in the two molecules.
    """
    properties_diff = {}
    deleted_atoms = set(atom.GetIdx() for atom in mol1.GetAtoms()) - set(atom_map.keys())
    for atom1_idx, atom2_idx in atom_map.items():
        atom1 = mol1.GetAtomWithIdx(atom1_idx)
        atom2 = mol2.GetAtomWithIdx(atom2_idx)

        # Check if the difference is related to intended changes
        if "added_hydrogens" in intended_changes:
            diff_implicit_valence = atom2.GetImplicitValence() - atom1.GetImplicitValence()
            diff_num_implicit_hs = atom2.GetNumImplicitHs() - atom1.GetNumImplicitHs()

            # If the difference in implicit valence and num implicit Hs is the same, ignore the difference
            if diff_implicit_valence == diff_num_implicit_hs:
                continue

        atom_diff = compare_atom_properties(atom1, atom2)
        # check if the difference is related to changed residue label, then remove from difference
        if "GetResidueName" in atom_diff and "changed_residue_label" in intended_changes:
            atom_diff.pop("GetResidueName")

        if atom_diff:
            properties_diff[atom1_idx] = atom_diff

    if "ignore_deleted" not in intended_changes:
        # Add deleted atoms to the differences
        for atom_idx in deleted_atoms:
            properties_diff[atom_idx] = {"DeletedAtom": True}

    return properties_diff


def get_properties(atom: Atom) -> dict:
    """Get properties of an RDKit atom."""
    res_dict = {}
    res = atom.GetPDBResidueInfo()
    res_dict['GetResidueName'] = res.GetResidueName()
    res_dict['GetResidueNumber'] = res.GetResidueNumber()
    res_dict['GetChainId'] = res.GetChainId()
    res_dict['GetName'] = res.GetName()
    prop_dict = {
        "GetSymbol": atom.GetSymbol(),
        "GetIdx": atom.GetIdx(),
        "GetAtomicNum": atom.GetAtomicNum(),
        "GetBonds": sorted([(sorted([bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()])) for bond in atom.GetBonds()]),
        "GetChiralTag": atom.GetChiralTag(),
        "GetDegree": atom.GetDegree(),
        "GetExplicitValence": atom.GetExplicitValence(),
        "GetFormalCharge": atom.GetFormalCharge(),
        "GetHybridization": atom.GetHybridization(),
        "GetImplicitValence": atom.GetImplicitValence(),
        "GetIsAromatic": atom.GetIsAromatic(),
        "GetIsotope": atom.GetIsotope(),
        "GetMass": atom.GetMass(),
        "GetNumExplicitHs": atom.GetNumExplicitHs(),
        "GetNumImplicitHs": atom.GetNumImplicitHs(),
        "IsInRing": atom.IsInRing()
    }
    prop_dict.update(res_dict)
    return prop_dict


def main() -> None:
    """Main function."""
    parser = argparse.ArgumentParser(description='Check if the topology has changed between two PDB files.')
    parser.add_argument('--file1', type=str, help='Path to the first PDB file')
    parser.add_argument('--file2', type=str, help='Path to the second PDB file')
    parser.add_argument('--intended_changes', nargs='+', type=str, default=[],
                        help='List of intended changes', choices=['added_hydrogens', 'changed_residue_label', 'ignore_deleted', 'atom_order_change'])
    args = parser.parse_args()
    mol1 = Chem.MolFromPDBFile(args.file1, removeHs=False, sanitize=False)
    mol2 = Chem.MolFromPDBFile(args.file2, removeHs=False, sanitize=False)
    signal.alarm(60)
    mol1_index_to_mol2_index = map_atoms(mol1, mol2)
    # Cancel the alarm if the code executed within the timeout duration
    signal.alarm(0)

    atom_differences = compare_atoms(mol1, mol2, mol1_index_to_mol2_index, intended_changes=args.intended_changes)
    if atom_differences:
        for mol1_idx, diff in atom_differences.items():
            if mol1_idx in mol1_index_to_mol2_index:  # otherwise was deleted
                mol2_idx = mol1_index_to_mol2_index[mol1_idx]
                print(f"File 1 Atom Idx {mol1_idx}, File 2 Atom Idx {mol2_idx}")
                for key, value in diff.items():
                    file1_value, file2_value = value
                    print(f"Property: {key} File 1: {file1_value}, File 2: {file2_value}")
            else:
                print(f"File 1 Atom Idx {mol1_idx} differences: {diff}")
        topology_change = True
    else:
        topology_change = False

    # NOTE: MDAnalysis changes original atom label and also atom order!!
    if 'atom_order_change' not in args.intended_changes:
        # Check if the atom order has changed
        if list(mol1_index_to_mol2_index.keys()) != list(mol1_index_to_mol2_index.values()):
            topology_change = True

    # Write the boolean value to the output file
    with open('valid.txt', 'w', encoding='utf-8') as file:
        file.write(str(topology_change))


if __name__ == '__main__':
    main()
