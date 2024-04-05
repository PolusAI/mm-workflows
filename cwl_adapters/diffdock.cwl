#!/usr/bin/env cwl-runner
cwlVersion: v1.0

class: CommandLineTool

label: DiffDock Diffusion based protein ligand docking

doc: |-
  DiffDock Diffusion based protein ligand docking

baseCommand: ["bash", "/DiffDock/diffdock_cmds.sh"]

hints:
  DockerRequirement:
    dockerPull: mrbrandonwalker/diffdock_cpu

requirements:
  InlineJavascriptRequirement: {}
inputs:

  protein_path:
    type: File
    format:
      - edam:format_1476
    inputBinding:
      prefix: --protein_path

  ligand_path:
    type: File
    format:
      - edam:format_3814
    inputBinding:
      prefix: --ligand_description

  inference_steps:
    label: number of reverse diffusion steps
    type: int?
    inputBinding:
      prefix: --inference_steps
    default: 20

  samples_per_complex:
    label: Number of sample poses to generate per complex
    type: int?
    inputBinding:
      prefix: --samples_per_complex
    default: 40

  batch_size:
    label: input batch size for neural net
    type: int?
    inputBinding:
      prefix: --batch_size
    default: 10

  out_dir:
    label: where output from diffdock is saved
    type: string?
    inputBinding:
      prefix: --out_dir
    default: results/

  model_dir:
    label: directory of DiffDock score model from paper
    type: string?
    inputBinding:
      prefix: --model_dir
    default:  /DiffDock/workdir/paper_score_model/

  confidence_model_dir:
    label: directory of DiffDock confidence model from paper
    type: string?
    inputBinding:
      prefix: --confidence_model_dir
    default:  /DiffDock/workdir/paper_confidence_model

  complex_name:
    label: name of folder with pose outputs that will be saved under out_dir folder
    type: string?
    inputBinding:
      prefix: --complex_name
    default:  outputs

outputs:

  max_confident_pose:
    type: File
    outputBinding:
      # the diffdock developers copy only the top ranked pose to a new file rank1.sdf
      glob: $(inputs.out_dir)/$(inputs.complex_name)/rank1.sdf
    format: edam:format_3814

  output_files:
    type: File[]
    outputBinding:
      # all other output files besides rank1.sdf have confidence information in them rank*_confidence*.sdf
      glob: $(inputs.out_dir)/$(inputs.complex_name)/rank*_confidence*.sdf
    format: edam:format_3814

  stderr:
    type: File
    outputBinding:
      glob: stderr

  execution_time:
    label: Time to run DiffDock
    doc: |-
      Time to run DiffDock
    type: float
    outputBinding:
      glob: stderr
      loadContents: true
      outputEval: |
        ${
          // the time command outputs to stderr and not to stdout
          // example output below, parse the float value of seconds (first item in line)
          // 1it [00:41, 41.03s/it]
          // 78.909
          return self[0].contents.split("\n").map(str => parseFloat(str)).reverse().find(num => !isNaN(num));
        }

stderr: stderr

$namespaces:
  edam: https://edamontology.org/
  cwltool: http://commonwl.org/cwltool#

$schemas:
- https://raw.githubusercontent.com/edamontology/edamontology/master/EDAM_dev.owl