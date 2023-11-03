#!/usr/bin/env cwl-runner
cwlVersion: v1.0

class: CommandLineTool

label: A tool that fixes structural issues in proteins

doc: |-
  A tool that fixes structural issues in proteins

baseCommand: ['python', '/pdb_fixer.py']

hints:
  DockerRequirement:
    dockerPull: ndonyapour/pdbfixer

inputs:
  input_pdb_path:
    type: File?
    format:
    - edam:format_1476
    inputBinding:
      prefix: --input_pdb_path

  input_helper_pdb_path:
    type: File?
    format:
    - edam:format_1476
    inputBinding:
      prefix: --input_helper_pdb_path

  pdbid:
    label: PDB ID from RCSB
    doc: |-
      PDB id from RCSB
      Type: string?
    type: string?
    format:
    - edam:format_2330 # 'Textual format'
    inputBinding:
      prefix: --pdbid

  url:
    label: URL to retrieve PDB from
    doc: |-
      URL to retrieve PDB from
      Type: string?
    type: string?
    format:
    - edam:format_2330 # 'Textual format'
    inputBinding:
      prefix: --url

  output_pdb_path:
    type: string
    format:
    - edam:format_1476
    inputBinding:
      prefix: --output_pdb_path
 
  add_atoms:
    label: What missing atoms to add, all, heavy or none
    doc: |-
      What missing atoms to add, all, heavy or none
      Type: string?
    type: string?
    format:
    - edam:format_2330 # 'Textual format'
    inputBinding:
      prefix: --add_atoms
    default: all

  add_residues:
    label: If set to True, adds missing residue
    doc: If set to True, adds missing residue
    type: boolean?
    format:
    - edam:format_2330
    inputBinding:
      prefix: --add_residues
    default: False

  replace_nonstandard:
    label: Replace nonstandard residues with standard equivalents
    doc: Replace nonstandard residues with standard equivalents
    type: boolean?
    format:
    - edam:format_2330
    inputBinding:
      prefix: --replace_nonstandard
    default: False

  keep_heterogens:
    # Note: A heterogen is any residue other than a standard amino acid or nucleotide
    label: What heterogens to keep, all, water or none 
    doc: What heterogens to keep, all, water or none 
    type: string?
    format:
    - edam:format_2330
    inputBinding:
      prefix: --keep_heterogens
    default: all

outputs:
  output_pdb_path:
    label: Output protein file path
    doc: |-
      Output protein file path
    type: File
    outputBinding:
      glob: $(inputs.output_pdb_path)
    format: edam:format_1476

$namespaces:
  edam: https://edamontology.org/

$schemas:
- https://raw.githubusercontent.com/edamontology/edamontology/master/EDAM_dev.owl
