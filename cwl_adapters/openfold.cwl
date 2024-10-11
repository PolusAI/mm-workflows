#!/usr/bin/env cwl-runner
cwlVersion: v1.2

class: CommandLineTool

label: OpenFold Inference

doc: |-
  OpenFold Inference using a pre-trained model.

baseCommand: ['python3', '/opt/openfold/run_pretrained_openfold.py']

hints:
  cwltool:CUDARequirement:
    cudaVersionMin: "10.2"
    cudaComputeCapability: "3.0"
    cudaDeviceCountMin: 1
    cudaDeviceCountMax: 1
  DockerRequirement:
    dockerPull: polusai/openfold-tool@sha256:84e1bb34f9d7ccfbd670e76d3f602334960f65ac091f96e9956ceeef985a9ffe
    dockerOutputDirectory: /app
  
requirements:
  EnvVarRequirement:
    envDef:
      CONDA_PREFIX: /opt/conda/ 
      # needed because https://github.com/aqlaboratory/openfold/blob/6f63267114435f94ac0604b6d89e82ef45d94484/scripts/utils.py#L9
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - $(inputs.input_fasta_dir)
      - $(inputs.template_mmcif_dir)
      - $(inputs.BASE_DATA_DIR)
      - $(inputs.PARAMS_DIR)

inputs:

  PARAMS_DIR:
    type: Directory
    doc: "Directory containing OpenFold model parameter files."

  BASE_DATA_DIR:
    type: Directory
    doc: "Directory containing OpenFold databases."

  input_fasta_dir:
    type: Directory
    doc: "Directory containing OpenFold input fasta sequence."

  template_mmcif_dir:
    type: Directory
    doc: "MMCIF files to use for template matching. This directory is required even if using template free inference."

  input_fasta_dir_string:
    type: string?
    inputBinding:
      position: 1
    default: /app/fasta

  template_mmcif_dir_string:
    type: string?
    inputBinding:
      position: 2
    default: /app/mmcif_files

  output_dir:
    type: string
    inputBinding:
      prefix: --output_dir
    default: "/app"
    doc: "Need to stage outputs in container in a specified directory"

  openfold_checkpoint_path:
    type: string?
    inputBinding:
      prefix: --openfold_checkpoint_path
    default: /app/params/finetuning_ptm_1.pt
    doc: "Uses an checkpoint or parameter file. Expected types are Deepspeed checkpoint files or .pt files. Make sure your selected checkpoint file matches the configuration setting chosen in --config_preset"

  config_preset:
    type: string?
    inputBinding:
      prefix: --config_preset
    default: "finetuning_ptm"
    doc: "Specify a different model configuration. There are 5 available model preset settings, some of which support template modeling, others support template-free modeling. The default is model_1."

  uniref90_database_path:
    type: string?
    inputBinding:
      prefix: --uniref90_database_path
    default: /app/openfold_databases/uniref90/uniref90.fasta
    doc: "For alignments to Uniref90 dataset."

  mgnify_database_path:
    type: string?
    inputBinding:
      prefix: --mgnify_database_path
    default: /app/openfold_databases/mgnify/mgy_clusters_2022_05.fa
    doc: "For alignments to  Mgnify, dataset."

  pdb70_database_path:
    type: string?
    inputBinding:
      prefix: --pdb70_database_path
    default: /app/openfold_databases/pdb70/pdb70
    doc: "For template matching to PDB70 database"

  uniclust30_database_path:
    type: string?
    inputBinding:
      prefix: --uniclust30_database_path
    default: /app/openfold_databases/uniclust30/uniclust30_2018_08/uniclust30_2018_08
    doc: "For template matching to Uniclust30 database"

  bfd_database_path:
    type: string?
    inputBinding:
      prefix: --bfd_database_path
    default: /app/openfold_databases/bfd/bfd_metaclust_clu_complete_id30_c90_final_seq.sorted_opt
    doc: "For alignments to BFD dataset."

  model_device:
    type: string?
    inputBinding:
      prefix: --model_device
    default: "cuda:0"

outputs:
  output_files:
    type: File[]
    outputBinding:
      glob: "/app/predictions/*"


$namespaces:
  edam: https://edamontology.org/
  cwltool: http://commonwl.org/cwltool#

$schemas:
- https://raw.githubusercontent.com/edamontology/edamontology/master/EDAM_dev.owl
