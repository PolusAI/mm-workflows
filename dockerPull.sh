#!/bin/bash -e
# Unlike the default github action runners which spin up a brand new machine every time,
# self-hosted runners are not necessarily isolated. In particular, the docker cache
# needs to be manually updated by explicitly performing docker pull commands.
# NOTE: For now use the following explicit list. In the future, consider using the
# cwl-utils library to recursively search for the dockerPull: tags within each workflow.
docker pull jakefennick/data
docker pull jakefennick/biosimspace

docker pull jakefennick/align_protein_ca_mda
docker pull jakefennick/align_protein_ca_pymol
docker pull jakefennick/atomselect
docker pull jakefennick/autodock_vina
docker pull jakefennick/autodock_vina_filter
docker pull jakefennick/bash_scripts
docker pull jakefennick/calculate_net_charge
docker pull jakefennick/generate_conformers
docker pull jakefennick/mol2_to_pdbqt
docker pull jakefennick/nmr4md
docker pull ndonyapour/openbabel
docker pull jakefennick/remove_terminal_residue_name_prefixes
docker pull jakefennick/rename_residues_mol
docker pull ndonyapour/molgan
docker pull cyangnyu/onionnet-sfct
docker pull cyangnyu/smina
docker pull ndonyapour/combine_structure
docker pull mrbrandonwalker/diffdock_gpu
docker pull mrbrandonwalker/diffdock_cpu
docker pull ndonyapour/pdbfixer
docker pull ndonyapour/extract_protein
docker pull mrbrandonwalker/extract_ligand_protein
docker pull ndonyapour/fix_pdb_atom_column
docker pull ndonyapour/generate_conformers