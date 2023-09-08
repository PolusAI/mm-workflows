#!/usr/bin/env cwl-runner
cwlVersion: v1.0

class: CommandLineTool

label: A tool that employs MDtraj to combine two PDB structures in a single PDB file.

doc: |-
  None

baseCommand: ['python', '/combine_pdb.py']

hints:
  DockerRequirement:
    dockerPull: ndonyapour/combine_pdb

inputs:
  input_structure1:
    label: Input structure 1 file path
    doc: |-
      Input structure 1 file path
      Type: string
      File type: input
      Accepted formats: pdb
      Example file: https://github.com/bioexcel/biobb_structure_utils/raw/master/biobb_structure_utils/test/data/utils/cat_protein.pdb
    type: File
    format:
    - edam:format_1476
    inputBinding:
      prefix: --input_structure1

  input_structure2:
    label: Input structure 2 file path
    doc: |-
      Input structure 2 file path
      Type: string
      File type: input
      Accepted formats: pdb
      Example file: https://github.com/bioexcel/biobb_structure_utils/raw/master/biobb_structure_utils/test/data/utils/cat_ligand.pdb
    type: File
    format:
    - edam:format_1476
    inputBinding:
      prefix: --input_structure2

  output_structure_path:
    label: Output protein file path
    doc: |-
      Output protein file path
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

  method:
    label: The frequency to save the outputs
    doc: |-
      The frequency to save the outputs
      Type: string
    type: string?
    format:
    - edam:format_2330
    inputBinding:
      prefix: --method
    default: mdtraj


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
