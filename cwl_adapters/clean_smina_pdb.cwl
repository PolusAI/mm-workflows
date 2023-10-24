#!/usr/bin/env cwl-runner
cwlVersion: v1.0

class: CommandLineTool

label: Clean smina pdb file (clean the pdb and rename the resname to LIG)

doc: |
  Clean smina pdb file (clean the pdb and rename the resname to LIG)

baseCommand: ["python", "/clean_smina_pdb.py"]

hints:
  DockerRequirement:
    dockerPull: cyangnyu/clean_smina_pdb

requirements:
  InlineJavascriptRequirement: {}

inputs:
  input_pdb:
    label: Input pdb file
    type: File
    format:
    - edam:format_1476
    inputBinding:
      prefix: --input_pdb

  output_pdb:
    label: Output pdb file
    type: string?
    format:
    - edam:format_1476
    inputBinding:
      prefix: --output_pdb

outputs:
  output_pdb:
    type: File
    format: edam:format_1476
    outputBinding:
      glob: $(inputs.output_pdb)

$namespaces:
  edam: https://edamontology.org/

$schemas:
- https://raw.githubusercontent.com/edamontology/edamontology/master/EDAM_dev.owl