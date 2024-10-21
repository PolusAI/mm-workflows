#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

label: OnionNet-SFCT tool for rescoring of docking poses

baseCommand: ["python", "/OnionNet-SFCT/scorer.py"]

hints:
  DockerRequirement:
    dockerPull: polusai/onionnet-sfct-tool@sha256:8b51dd5e9d1218a4033cb6126c44aed0ea67a7ad4cfcf49ee4703fce7520c6be

requirements:
  InlineJavascriptRequirement: {}

inputs:

  receptor_path:
    label: receptor structure in in pdb, pdbqt, mol, mol2, sdf format
    type: File?
    format:
    - edam:format_1476
    - edam:format_3814
    - edam:format_3815
    - edam:format_3816
    inputBinding:
      prefix: -r

  ligand_path:
    label: ligand docking poses structure in pdb or pdbqt format.
    type: File?
    format: edam:format_1476
    inputBinding:
      prefix: -l

  model_path:
    label: OnionNet-SFCT final model
    type: string?
    inputBinding:
      prefix: --model
    default: "/OnionNet-SFCT/data/sfct_std_final.model"

  pose_type:
    label: docking poses type, can be ["vina", "smina", "gnina", "idock", "ledock"]
    type: string?
    inputBinding:
      prefix: --stype
    default: "vina"

  output_score_path:
    label: text file contains docking_score and onionnet_rescore for each pose
    type: string?
    inputBinding:
      prefix: -o
    default: "sfct.txt"

# Is this still necessary?
  output_path:
    type: string

outputs:

  output_score_path:
    type: File
    outputBinding:
      glob: $(inputs.output_score_path)

  output_path:
    type: File
    outputBinding:
      glob: $(inputs.output_path)

  output_docking_score:
    label: Estimated Free Energy of Binding (docking score)
    doc: |-
      Estimated Free Energy of Binding
    type: float
    outputBinding:
      glob: $(inputs.output_score_path)
      loadContents: true
      outputEval: |
        ${
          const lines = self[0].contents.split("\n");
          // The correct line should be of the form
          // # name pose_index origin_score combined_score sfct
          // docked_vina.pdb 0 -10.365 -4.137 2.090
          // pose_1 1 -9.973 -3.509 2.956
          const bfe_line = lines[1];
          const docking_score_string = bfe_line.split(" ").filter(function(s) {return !isNaN(parseFloat(s))})[1];
          const output_docking_score = parseFloat(docking_score_string);
          return output_docking_score
        }

  output_poses_rescore:
    label: Estimated Free Energy of Binding (onionnet rescore)
    doc: |-
      Estimated Free Energy of Binding
    type: float
    outputBinding:
      glob: $(inputs.output_score_path)
      loadContents: true
      outputEval: |
        ${
          const lines = self[0].contents.split("\n");
          // The correct line should be of the form
          // # name pose_index origin_score combined_score sfct
          // docked_vina.pdb 0 -10.365 -4.137 2.090
          // pose_1 1 -9.973 -3.509 2.956
          const bfe_line = lines[1];
          const docking_score_string = bfe_line.split(" ").filter(function(s) {return !isNaN(parseFloat(s))})[2];
          const output_poses_rescore = parseFloat(docking_score_string);
          return output_poses_rescore
        }

stdout: $(inputs.output_path)

$namespaces:
  edam: https://edamontology.org/
  cwltool: http://commonwl.org/cwltool#

$schemas:
- https://raw.githubusercontent.com/edamontology/edamontology/master/EDAM_dev.owl
