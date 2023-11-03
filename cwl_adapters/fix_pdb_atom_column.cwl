#!/usr/bin/env cwl-runner
cwlVersion: v1.0

class: CommandLineTool

label: Modifies the atom element types using a helper structure.

doc: |-
  Modifies the atom element types using a helper structure.

baseCommand: ['python3', '/fix_pdb_atom_column.py']

hints:
  DockerRequirement:
    dockerPull: ndonyapour/fix_pdb_atom_column

inputs:

  input_structure_path:
    type: File
    format:
    - edam:format_1476
    inputBinding:
      position: 1

  input_helper_structure_path:
    type: File
    format:
    - edam:format_1476
    inputBinding:
      position: 2

  output_pdb_path:
    type: string
    format:
    - edam:format_1476
    inputBinding:
     position: 3
    default: system.pdb

outputs:
  output_pdb_path:
    type: File
    format: edam:format_1476
    streamable: true
    outputBinding:
      glob: $(inputs.output_pdb_path)

stdout: $(inputs.output_pdb_path)

$namespaces:
  edam: https://edamontology.org/

$schemas:
- https://raw.githubusercontent.com/edamontology/edamontology/master/EDAM_dev.owl
