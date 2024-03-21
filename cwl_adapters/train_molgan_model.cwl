#!/usr/bin/env cwl-runner
cwlVersion: v1.0

class: CommandLineTool

label: MolGAN tool for generating small molecules

baseCommand: ["python", "/MolGAN/example.py"]

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

  output_model_dir:
    label: Output directory
    type: string
    format:
    - edam:format_2330 # 'Textual format'
    inputBinding:
      prefix: --output_model_dir
    default: output

  validation_metrics:
    label: The metrics are used during validation and testing. Metrics, 'np,logp,sas,qed,novelty,dc,unique,diversity,validity'
    doc: |-
      The metrics are used during validation and testing
    type: string?
    format:
    - edam:format_2330
    inputBinding:
      prefix: --validation_metrics
    default: 'np,logp,sas,qed,novelty,dc,unique,diversity,validity'

  num_epochs:
    label: The number of training epochs
    doc: |-
      The number of training epochs
      Type: int
    type: int?
    format:
    - edam:format_2330
    inputBinding:
      prefix: --num_epochs

  save_frequency:
    label: The frequency to save the outputs
    doc: |-
      The frequency to save the outputs
      Type: int
    type: int?
    format:
    - edam:format_2330
    inputBinding:
      prefix: --save_frequency

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

  output_model_dir:
    type: Directory
    outputBinding:
      glob: $(inputs.output_model_dir)
    format: edam:format_2330 # 'Textual format

  stderr:
    type: File
    outputBinding:
      glob: stderr

stderr: stderr

$namespaces:
  edam: https://edamontology.org/

$schemas:
- https://raw.githubusercontent.com/edamontology/edamontology/master/EDAM_dev.owl

