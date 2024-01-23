#!/usr/bin/env cwl-runner
cwlVersion: v1.0

class: CommandLineTool

label: P2rank protein pocket detection and ranking

doc: |-
  P2rank protein pocket detection and ranking

baseCommand: ["/p2rank/prank", "predict"]

hints:
  DockerRequirement:
    dockerPull: mrbrandonwalker/p2rank

requirements:
  InlineJavascriptRequirement: {}

inputs:

  protein_filename:
    type: File
    format:
      - edam:format_1476
    inputBinding:
      prefix: -f

  out_dir:
    type: string
    inputBinding:
      prefix: -o
    default: "test_output"

  pocket_predictions:
    type: string?

outputs:

  pocket_predictions:
    type: File
    outputBinding:
      # If input filename is 1fbl.pdb, then output filename is 1fbl.pdb_predictions.csv
      glob: $(inputs.out_dir)/$(inputs.protein_filename.basename)_predictions.csv
    format: edam:format_3752

$namespaces:
  edam: https://edamontology.org/
  cwltool: http://commonwl.org/cwltool#

$schemas:
- https://raw.githubusercontent.com/edamontology/edamontology/master/EDAM_dev.owl