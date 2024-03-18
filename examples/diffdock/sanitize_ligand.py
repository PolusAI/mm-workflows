# mypy: disable-error-code="attr-defined"
# pylint: disable=E0401,E1101,I1101
"""Filter molecules with kekulization errors.
Handle molecules with sanitization errors by assign formal charge based on valence."""
# https://depth-first.com/articles/2020/02/10/a-comprehensive-treatment-of-aromaticity-in-the-smiles-language/
import argparse
from pathlib import Path
from typing import Dict, Tuple, Optional
from types import FrameType
import signal  # only unix supported

import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem

parser = argparse.ArgumentParser()
parser.add_argument('--input_small_mol_ligand', type=str, help='Ligand file')
args = parser.parse_args()
input_small_mol_ligand = Path(args.input_small_mol_ligand).name


class TimeoutException(Exception):
    pass


def timeout(seconds: int, error_message: str = 'Timeout') -> None:
    """Time out for conformer generation, only works for unix. For windows have to ignore mypy error.

    Args:
        seconds (int): The time in seconds to wait before timing out
        error_message (str, optional): The error message to display. Defaults to 'Timeout'.
    """
    def handle_timeout(signum: int, frame: Optional[FrameType]) -> None:
        raise TimeoutException(error_message)

    signal.signal(signal.SIGALRM, handle_timeout)
    signal.alarm(seconds)


def adjust_formal_charges(molecule: Chem.SDMolSupplier) -> Chem.SDMolSupplier:
    """Sometimes input structures do not have correct formal charges corresponding
    to bond order topology. So choose to trust bond orders assigned and generate formal
    charges based on that.
    Explicit valence determined what the formal charge should be from dictionary of valence
    to formal charge for that atomic number. Special case if atom == carbon or nitrogen
    and if neighbors contain nitrogen, oyxgen or sulfur (polarizable atoms) then if carbon
    and explicit valence only 3, give formal charge of +1 (more stable then -1 case).

    Args:
        molecule (Chem.SDMolSupplier): The rdkit molecule object

    Returns:
        Chem.SDMolSupplier: Molecule object with adjusted formal charges
    """
    # 1=H, 5=B, 6=C, 7=N, 8=O, 15=P, 16=S, 17=Cl, 9=F, 35=Br, 53=I
    atomicnumtoformalchg: Dict[int, Dict[int, int]] = {1: {2: 1}, 5: {4: 1}, 6: {3: -1}, 7: {2: -1, 4: 1},
                                                       8: {1: -1, 3: 1}, 15: {4: 1}, 16: {1: -1, 3: 1, 5: -1},
                                                       17: {0: -1, 4: 3}, 9: {0: -1}, 35: {0: -1}, 53: {0: -1}}
    for atom in molecule.GetAtoms():
        atomnum = atom.GetAtomicNum()
        val = atom.GetExplicitValence()
        if atomicnumtoformalchg.get(atomnum) is None:
            continue
        valtochg = atomicnumtoformalchg[atomnum]
        chg = valtochg.get(val, 0)
        # special case of polar neighbors surrounding carbon or nitrogen
        # https://docs.eyesopen.com/toolkits/cpp/oechemtk/valence.html#openeye-charge-model
        #
        # See Section 6: Factors That Stabilize Carbocations – Adjacent Lone Pairs
        # https://www.masterorganicchemistry.com/2011/03/11/3-factors-that-stabilize-carbocations/#six
        polneighb = False
        if atomnum in (6, 7):
            for natom in atom.GetNeighbors():
                natomicnum = natom.GetAtomicNum()
                if natomicnum in (7, 8, 16):
                    polneighb = True
            if polneighb and val == 3 and atomnum == 6:
                chg = 1

        atom.SetFormalCharge(chg)
    return molecule

# for some ligands can take up to minutes to generate conformer!!


def generate_conformer(molecule: Chem.SDMolSupplier, timeout_value: int = 5) -> None:
    """ Generate conformer for molecule. Sometimes rdkit embedding can fail,
      so use random coordinates if failed at first.

    Args:
        molecule (Chem.SDMolSupplier): The rdkit molecule object
        timeout_value (float): The time in seconds to wait before timing out
    """
    # Set up the signal handler for the timeout
    timeout(timeout_value)
    ps = AllChem.ETKDGv2()
    confid = AllChem.EmbedMolecule(molecule, ps)
    if confid == -1:
        ps.useRandomCoords = True
        AllChem.EmbedMolecule(molecule, ps)
        # only want to catch error Bad Conformation Id,  dont need to spend time optimizing
        AllChem.MMFFOptimizeMolecule(molecule, confId=0, maxIters=1)

    signal.alarm(0)


def is_valid_ligand(molecule: Chem.SDMolSupplier) -> Tuple[bool, Chem.SDMolSupplier]:
    """Check for sanitization errors, attempt to fix formal charges/valence consistency errors.
       DiffDock uses rdkit to generate a seed conformation that will sometimes crash, so generating
       conformations here to catch that error and prevent DiffDock from running that ligand.

    Args:
        mol (Chem.SDMolSupplier): The rdkit small molecule object

    Returns:
        bool: if ligand is valid
        Chem.SDMolSupplier: molecule object
    """
    valid_lig = True
    try:
        Chem.SanitizeMol(molecule)
        generate_conformer(molecule)
    except rdkit.Chem.rdchem.KekulizeException as e:
        valid_lig = False
        # Not handling kekulization error now so just remove file to prevent DiffDock execution
    except rdkit.Chem.rdchem.MolSanitizeException as e:
        # can also be explicit valence error (i.e.) formal charge not consistent with bond topology
        # choose to trust bond topology around atom and add formal charge based on that
        molecule = adjust_formal_charges(molecule)
    except ValueError:
        # assuming this is Bad Conformer Id error from generate_conformer
        valid_lig = False
    except Exception:  # pylint: disable=broad-exception-caught
        # catch *all* exceptions rdkit can throw
        valid_lig = False

    return valid_lig, molecule


mol: Chem.SDMolSupplier = Chem.SDMolSupplier(args.input_small_mol_ligand, sanitize=False, removeHs=False)[0]

valid_ligand, rdkit_mol = is_valid_ligand(mol)

if valid_ligand:
    with Chem.SDWriter(input_small_mol_ligand) as w:
        w.write(rdkit_mol)
