#!/usr/bin/env cwl-runner
cwlVersion: v1.0

class: CommandLineTool

label: A tool that employs OpenMM to extract protein from a PDB file

doc: |-
  A tool that employs OpenMM to extract protein from a PDB file

baseCommand: ['python', '/extract_protein.py']

hints:
  DockerRequirement:
    dockerPull: ndonyapour/extract_protein

inputs:
  input_pdb_path:
    label: Input pdb file path
    doc: |-
      Input pdb file path
      Type: string
      File type: input
      Accepted formats: pdb
      Example file: https://github.com/bioexcel/biobb_structure_utils/raw/master/biobb_structure_utils/test/data/utils/cat_protein.pdb
    type: File
    format:
    - edam:format_1476
    inputBinding:
      prefix: --input_pdb_path

  output_pdb_path:
    label: Output pdb file path
    doc: |-
      Output pdb file path
      Type: string
      File type: output
      Accepted formats: pdb
      Example file: https://github.com/bioexcel/biobb_structure_utils/raw/master/biobb_structure_utils/test/reference/utils/ref_cat_pdb.pdb
    type: string
    format:
    - edam:format_1476
    inputBinding:
      prefix: --output_pdb_path
    default: system.pdb

outputs:
  output_pdb_path:
    label: Output pdb file path
    doc: |-
      Output pdb file path
    type: File
    outputBinding:
      glob: $(inputs.output_pdb_path)
    format: edam:format_1476

$namespaces:
  edam: https://edamontology.org/

$schemas:
- https://raw.githubusercontent.com/edamontology/edamontology/master/EDAM_dev.owl