#!/usr/bin/env cwl-runner
cwlVersion: v1.0

class: CommandLineTool

label: Change the labels of residues

doc: |
  Change the labels of residues

baseCommand: ["python", "/rename_residues.py"]

hints:
  DockerRequirement:
    dockerPull: ndonyapour/rename_residues

inputs:
  input_path:
    type: File
    format:
    - edam:format_3816 # mol2
    inputBinding:
      prefix: --input_path

  output_path:
    label: Path to the output file
    doc: |-
      Path to the output file
    type: string
    format:
    - edam:format_3816 # mol2
    inputBinding:
      prefix: --output_path
    default: system.mol2

  residue_name:
    label: The new residue name to which all residues should be changed 
    doc: |-
      The new residue name to which all residues should be changed
    type: string
    format:
    - edam:format_2330
    inputBinding:
      prefix: --residue_name 
 
  column_index:
    label: The index of column of the line
    doc: |-
     The index of column of the line
    type: int?
    format:
    - edam:format_2330
    inputBinding:
      prefix: --column_index

  line_start_column_idxs:
    label: A dictionary string indicating the line start and corresponding column indice 
    doc: |-
      A dictionary string indicating the line start and corresponding column indice
    type: string?
    format:
    - edam:format_2330
    inputBinding:
      prefix: --line_start_column_idxs
outputs:
  output_path:
    label: Path to the output file
    doc: |-
      Path to the output file
    type: File
    format: edam:format_3816 # mol2
    outputBinding:
      glob: $(inputs.output_path)

$namespaces:
  edam: https://edamontology.org/

$schemas:
- https://raw.githubusercontent.com/edamontology/edamontology/master/EDAM_dev.owl