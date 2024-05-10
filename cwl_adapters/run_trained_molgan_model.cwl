#!/usr/bin/env cwl-runner
cwlVersion: v1.0

class: CommandLineTool

label: MolGAN tool for generating small molecules

baseCommand: ["python", "/MolGAN/run_trained_model.py"]

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
  input_data_path:
    label: Path to the input data file
    doc: |-
      Path to the input data file
      Type: File
      File type: input
      Accepted formats: pkl
      Example file: https://github.com/bioexcel/biobb_ml/raw/master/biobb_ml/test/reference/classification/ref_output_model_support_vector_machine.pkl
    type: File
    format:
    - edam:format_3653
    inputBinding:
      prefix: --input_data_path
    default: system.pkl

  input_NP_Score_path:
    label: Output ceout file (AMBER ceout)
    doc: |-
      Output ceout file (AMBER ceout)
      Type: File
      File type: input
      Accepted formats: gz
      Example file: https://github.com/bioexcel/biobb_amber/raw/master/biobb_amber/test/data/cphstats/sander.ceout.gz
    type: File
    format:
    - edam:format_3987
    default: NP.gz
    inputBinding:
      prefix: --input_NP_Score_path

  input_SA_Score_path:
    label: Output ceout file (AMBER ceout)
    doc: |-
      Output ceout file (AMBER ceout)
      Type: File
      File type: input
      Accepted formats: gz
      Example file: https://github.com/bioexcel/biobb_amber/raw/master/biobb_amber/test/data/cphstats/sander.ceout.gz
    type: File
    format:
    - edam:format_3987
    default: SA.gz
    inputBinding:
      prefix: --input_SA_Score_path

  input_model_dir:
    label: Input directory of trained models
    type: Directory
    format:
    - edam:format_2330 # 'Textual format'
    inputBinding:
      prefix: --input_model_dir
    default: output

  output_log_path:
    label: Path to the log file
    doc: |-
      Path to the log file
      Type: string
      File type: output
      Accepted formats: log
    type: string
    format:
    - edam:format_2330
    inputBinding:
      prefix: --output_log_path
    default: system.log

  output_sdf_path:
    label: Path to the output file
    doc: |-
      Path to the output file
      Type: string
      File type: output
      Accepted formats: sdf
    type: string
    format:
    - edam:format_3814 # sdf
    default: system.sdf
    inputBinding:
      prefix: --output_sdf_path

  num_samples:
    label: The number of new molecules to generate
    doc: |-
      The number of training epochs
      Type: int
    type: int?
    format:
    - edam:format_2330
    inputBinding:
      prefix: --num_samples
    default: 1000

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
  output_log_path:
    label: Path to the log file
    doc: |-
      Path to the log file
    type: File
    outputBinding:
      glob: $(inputs.output_log_path)
    format: edam:format_2330

  output_sdf_path:
    label: Path to the output file
    doc: |-
      Path to the output file
    type: File
    outputBinding:
      glob: $(inputs.output_sdf_path)
    format: edam:format_3814 # sdf

  stderr:
    type: File
    outputBinding:
      glob: stderr

stderr: stderr


$namespaces:
  edam: https://edamontology.org/

$schemas:
- https://raw.githubusercontent.com/edamontology/edamontology/master/EDAM_dev.owl

