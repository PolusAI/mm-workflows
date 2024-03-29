#!/usr/bin/env cwl-runner
cwlVersion: v1.0

class: CommandLineTool

label: Check for topology changes in PDB files

doc: |-
  Check for topology changes in PDB files

baseCommand: ["python", "/topology_check.py"]

hints:
  DockerRequirement:
    dockerPull: mrbrandonwalker/topology_check

requirements:
  InlineJavascriptRequirement: {}

inputs:

  file1:
    type: File
    inputBinding:
      prefix: --file1
    format:
    - edam:format_1476

  file2:
    type: File
    inputBinding:
      prefix: --file2
    format:
    - edam:format_1476

  intended_changes:
    type: string?
    inputBinding:
      prefix: --intended_changes

  topology_changed:
    type: string?

outputs:

  topology_changed:
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

  stdout:
    type: File
    outputBinding:
      glob: stdout

stdout: stdout

$namespaces:
  edam: https://edamontology.org/

$schemas:
- https://raw.githubusercontent.com/edamontology/edamontology/master/EDAM_dev.owl