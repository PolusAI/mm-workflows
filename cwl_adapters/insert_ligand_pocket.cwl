#!/usr/bin/env cwl-runner
cwlVersion: v1.0

class: CommandLineTool

label: Insert ligand into protein pocket(s)

doc: |-
  Insert ligand into protein pocket(s)

baseCommand: ["python", "/insert_ligand_pocket.py"]

hints:
  DockerRequirement:
    dockerPull: mrbrandonwalker/insert_ligand_pocket

requirements:
  InlineJavascriptRequirement: {}

inputs:

  protein_filename:
    type: File
    format:
      - edam:format_1476
    inputBinding:
      prefix: --protein_filename

  ligand_filename:
    type: File
    format:
      - edam:format_3814
    inputBinding:
      prefix: --ligand_filename

  input_pocket_predictions:
    type: File
    format: edam:format_3752
    inputBinding:
      prefix: --pocket_predictions

  top_n_pockets:
    type: int
    inputBinding:
      prefix: --top_n_pockets

  protein_ligand_complexes:
    type: string?

outputs:

  protein_ligand_complexes:
    type: File[]
    outputBinding:
      glob: "*.pdb"
    format: edam:format_1476

$namespaces:
  edam: https://edamontology.org/
  cwltool: http://commonwl.org/cwltool#

$schemas:
- https://raw.githubusercontent.com/edamontology/edamontology/master/EDAM_dev.owl