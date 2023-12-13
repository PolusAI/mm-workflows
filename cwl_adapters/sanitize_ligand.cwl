#!/usr/bin/env cwl-runner
cwlVersion: v1.0

class: CommandLineTool

label: Sanitize input ligand

doc: |-
  Sanitize input ligand

baseCommand: ["python", "/sanitize_ligand.py"]

hints:
  DockerRequirement:
    dockerPull:  mrbrandonwalker/sanitize_ligand

requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement: # conditionally overwrite the input ligand, otherwise cwltool will symlink to the original
    listing:
      - $(inputs.input_small_mol_ligand)

inputs:

  input_small_mol_ligand:
    type: File
    format:
      - edam:format_3814
    inputBinding:
      prefix: --input_small_mol_ligand

  output_ligand:
    type: string?

  valid_ligand:
    type: string?

outputs:

  output_ligand:
    type: File
    format: edam:format_3814
    outputBinding:
      glob: "*.sdf"

  valid_ligand:
    type: boolean
    outputBinding:
      glob: valid.txt
      loadContents: true
      outputEval: |
        ${
          // Read the contents of the file
          const lines = self[0].contents.split("\n");
          // Read boolean value from the first line
          const valid = lines[0].trim() === "True";
          return valid;

        }

$namespaces:
  edam: https://edamontology.org/

$schemas:
- https://raw.githubusercontent.com/edamontology/edamontology/master/EDAM_dev.owl