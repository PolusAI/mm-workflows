#!/usr/bin/env cwl-runner
cwlVersion: v1.0

class: CommandLineTool

label: A tool that employs MDAnalysis to combine two PDB structures in a single PDB file.

doc: |-
  A tool that employs MDAnalysis to combine two PDB structures in a single PDB file.

baseCommand: ['python', '/combine_pdb.py']

hints:
  DockerRequirement:
    dockerPull: ndonyapour/combine_pdb

inputs:
  input_structure1:
    label: Input PDB structure 1 file path
    doc: |-
      Input PDB structure 1 file path
      Type: string
      File type: input
      Accepted formats: pdb
    type: File
    format:
    - edam:format_1476
    inputBinding:
      prefix: --input_structure1

  input_structure2:
    label: Input PDB structure 2 file path
    doc: |-
      Input PDB structure 2 file path
      Type: string
      File type: input
      Accepted formats: pdb
    type: File
    format:
    - edam:format_1476
    inputBinding:
      prefix: --input_structure2

  output_structure_path:
    label: Output combined PDB file path
    doc: |-
      Output combined PDB file path
      Type: string
      File type: output
      Accepted formats: pdb
      Example file: https://github.com/bioexcel/biobb_structure_utils/raw/master/biobb_structure_utils/test/reference/utils/ref_cat_pdb.pdb
    type: string
    format:
    - edam:format_1476
    inputBinding:
      prefix: --output_structure_path
    default: system.pdb
outputs:
  output_structure_path:
    label: Output protein file path
    doc: |-
      Output protein file path
    type: File
    outputBinding:
      glob: $(inputs.output_structure_path)
    format: edam:format_1476

$namespaces:
  edam: https://edamontology.org/

$schemas:
- https://raw.githubusercontent.com/edamontology/edamontology/master/EDAM_dev.owl
