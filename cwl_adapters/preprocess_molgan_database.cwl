#!/usr/bin/env cwl-runner
cwlVersion: v1.0

class: CommandLineTool

label: MolGAN tool for generating small molecules

baseCommand: ["python", "/MolGAN/utils/sparse_molecular_dataset.py"]

hints:
  DockerRequirement:
    dockerPull: ndonyapour/molgan

# Set environment variables for the tool,
# See: https://www.commonwl.org/user_guide/topics/environment-variables.html
requirements:
  EnvVarRequirement:
    envDef:
      RDKIT_ERROR_LOGGING: $(inputs.rdkit_error_logging)

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

  rdkit_error_logging:
    label: Enable or disable RDKit error logging
    doc: |-
      Enable or disable RDKit error logging
    type: string?
    format:
    - edam:format_2330
    # RDKit prints out all errors by default, which can pose issues for CI,
    # particularly with large databases. It would be more efficient to suppress these errors.
    default: "ON"

outputs:
  output_data_path:
    label: Path to the output data file
    doc: |-
      Path to the output data file
    type: File
    outputBinding:
      glob: $(inputs.output_data_path)
    format: edam:format_3653 # sdf

  stderr:
    type: File
    outputBinding:
      glob: stderr

stderr: stderr

$namespaces:
  edam: https://edamontology.org/

$schemas:
- https://raw.githubusercontent.com/edamontology/edamontology/master/EDAM_dev.owl

