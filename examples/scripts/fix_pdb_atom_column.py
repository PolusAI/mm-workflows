import sys

import MDAnalysis as mda  # pylint: disable=import-error

input_structure = sys.argv[1]
input_helper_structure = sys.argv[2]
output_pdb = sys.argv[3]

u1 = mda.Universe(input_structure)
u2 = mda.Universe(input_helper_structure)

# replace the atom elements
u1.atoms.elements = u2.atoms.elements

all_atoms = u1.select_atoms("all")
all_atoms.write(output_pdb)
