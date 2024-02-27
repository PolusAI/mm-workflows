#!/usr/bin/env cwl-runner
cwlVersion: v1.0

class: CommandLineTool

label: MolGAN tool for generating small molecules

baseCommand: ["python", "/MolGAN/utils/sparse_molecular_dataset.py"]

hints:
  DockerRequirement:
    dockerPull: ndonyapour/molgan

inputs:
  input_sdf_path:
    label: Path to the input file
    doc: |-
      Path to the input file
      Type: File
      File type: input
      Accepted formats: sdf
    type: File
    format:
    - edam:format_3814 # sdf
    default: system.sdf
    inputBinding:
      prefix: --input_sdf_path

  output_data_path:
    label: Path to the output data file
    doc: |-
      Path to the output data file
      Type: string
      File type: output
      Accepted formats: pkl
      Example file: https://github.com/bioexcel/biobb_ml/raw/master/biobb_ml/test/reference/classification/ref_output_model_support_vector_machine.pkl
    type: string
    format:
    - edam:format_3653
    inputBinding:
      prefix: --output_data_path
    default: system.pkl

outputs:
  output_data_path:
    label: Path to the output data file
    doc: |-
      Path to the output data file
    type: File
    outputBinding:
      glob: $(inputs.output_data_path)
    format: edam:format_3653 # sdf

$namespaces:
  edam: https://edamontology.org/

$schemas:
- https://raw.githubusercontent.com/edamontology/edamontology/master/EDAM_dev.owl

