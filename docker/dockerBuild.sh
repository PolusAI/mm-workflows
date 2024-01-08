#!/bin/bash -e
# Unlike the default github action runners which spin up a brand new machine every time,
# self-hosted runners are not necessarily isolated. In particular, the docker cache
# needs to be manually updated by explicitly performing docker pull/build commands.
# NOTE: For now use the following explicit list. In the future, consider using the
# cwl-utils library to recursively search for the dockerPull: tags within each workflow.

cd ..

cd examples/data/
sudo docker build --no-cache --pull -f Dockerfile_data -t jakefennick/data .
cd ../..

cd cwl_adapters/file_format_conversions/biosimspace/
sudo docker build --no-cache --pull -f Dockerfile_biosimspace -t jakefennick/biosimspace .
cd ../../..

cd examples/scripts/
sudo docker build --no-cache --pull -f Dockerfile_align_protein_ca_mda -t jakefennick/align_protein_ca_mda .
sudo docker build --no-cache --pull -f Dockerfile_align_protein_ca_pymol -t jakefennick/align_protein_ca_pymol .
sudo docker build --no-cache --pull -f Dockerfile_atomselect -t jakefennick/atomselect .
sudo docker build --no-cache --pull -f Dockerfile_autodock_vina -t jakefennick/autodock_vina .
sudo docker build --no-cache --pull -f Dockerfile_autodock_vina_filter -t jakefennick/autodock_vina_filter .
sudo docker build --no-cache --pull -f Dockerfile_bash_scripts -t jakefennick/bash_scripts .
sudo docker build --no-cache --pull -f Dockerfile_calculate_net_charge -t jakefennick/calculate_net_charge .
sudo docker build --no-cache --pull -f Dockerfile_mol2_to_pdbqt -t jakefennick/mol2_to_pdbqt .
sudo docker build --no-cache --pull -f Dockerfile_nmr4md -t jakefennick/nmr4md .
sudo docker build --no-cache --pull -f Dockerfile_openbabel -t jakefennick/openbabel .
sudo docker build --no-cache --pull -f Dockerfile_remove_terminal_residue_name_prefixes -t jakefennick/remove_terminal_residue_name_prefixes .
sudo docker build --no-cache --pull -f Dockerfile_rename_residues_mol -t jakefennick/rename_residues_mol .
sudo docker build --no-cache --pull -f Dockerfile_combine_structure -t ndonyapour/combine_structure .
sudo docker build --no-cache --pull -f Dockerfile_molgan -t ndonyapour/molgan .
sudo docker build --no-cache --pull -f Dockerfile_onionnet-sfct -t cyangnyu/onionnet-sfct .
sudo docker build --no-cache --pull -f Dockerfile_smina -t cyangnyu/smina .
sudo docker build --no-cache --pull -f Dockerfile_pdb_fixer -t ndonyapour/pdbfixer .
sudo docker build --no-cache --pull -f Dockerfile_extract_protein -t ndonyapour/extract_protein .
sudo docker build --no-cache --pull -f Dockerfile_fix_pdb_atom_column -t ndonyapour/fix_pdb_atom_column .
sudo docker build --no-cache --pull -f Dockerfile_generate_conformers -t ndonyapour/generate_conformers .

sudo docker build --no-cache --pull -f Dockerfile_pdbbind_refined -t pdbbind_refined_v2020 .  # NOTE: no username
cd ../..

cd examples/diffdock/
sudo docker build --no-cache --pull -f Dockerfile_diffdock_cpu -t mrbrandonwalker/diffdock_cpu .
sudo docker build --no-cache --pull -f Dockerfile_diffdock_gpu -t mrbrandonwalker/diffdock_gpu .
sudo docker build --no-cache --pull -f Dockerfile_rmsd_pose_cluster -t mrbrandonwalker/rmsd_pose_cluster .
sudo docker build --no-cache --pull -f Dockerfile_rank_diffdock_poses -t mrbrandonwalker/rank_diffdock_poses .
cd ../..