#!/usr/bin/env cwl-runner
cwlVersion: v1.0

class: CommandLineTool

label: Visualize output molecular structure files in pymol

doc: |-
  Visualize output molecular structure files in pymol

baseCommand: ["python", "/visualize_pymol.py"]

hints:
  DockerRequirement:
    dockerPull: mrbrandonwalker/visualize_pymol

requirements:
  InlineJavascriptRequirement: {}

inputs:

  pdb_path:
    type: File
    inputBinding:
      prefix: --pdb_path

  ligand_path:
    type: File?
    inputBinding:
      prefix: --ligand_path

  output_pymol_session:
    type: string?

outputs:

  output_pymol_session:
    type: File
    outputBinding:
      glob: "*.pse"

$namespaces:
  edam: https://edamontology.org/

$schemas:
- https://raw.githubusercontent.com/edamontology/edamontology/master/EDAM_dev.owl