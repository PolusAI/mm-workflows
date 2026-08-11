#!/usr/bin/env cwl-runner
cwlVersion: v1.0

class: CommandLineTool

label: Run a Bash script

doc: |
  Run a Bash script

# NOTE: The `ndonyapour/molgan` image's baked-in /MolGAN/data/download_dataset.sh
# downloads the gdb9 dataset from http://deepchem.io.s3-website-us-west-1.amazonaws.com,
# which now returns 403 Forbidden / AllAccessDisabled (DeepChem retired that old
# S3-website alias). This overrides that script with a corrected copy pointing at
# DeepChem's current canonical bucket (verified against deepchem/deepchem's own
# qm9_datasets.py loader: GDB9_URL = "https://deepchemdata.s3-us-west-1.amazonaws.com/...")
# instead of waiting on an upstream fix/rebuild of the third-party image.
baseCommand: ["bash", "-c"]

arguments:
- position: 1
  valueFrom: |
    set -euo pipefail
    wget -nv --no-clobber https://deepchemdata.s3-us-west-1.amazonaws.com/datasets/gdb9.tar.gz
    tar xvzf gdb9.tar.gz
    rm gdb9.tar.gz
    rm gdb9.sdf.csv
    if [ "gdb9.sdf" != "$1" ]; then
       mv gdb9.sdf "$1"
    fi
    wget -nv --no-clobber https://github.com/gablg1/ORGAN/raw/master/organ/NP_score.pkl.gz
    wget -nv --no-clobber https://github.com/gablg1/ORGAN/raw/master/organ/SA_score.pkl.gz
    if [ "NP_score.pkl.gz" != "$2" ]; then
       mv NP_score.pkl.gz "$2"
    fi
    if [ "SA_score.pkl.gz" != "$3" ]; then
       mv SA_score.pkl.gz "$3"
    fi
- position: 2
  valueFrom: download_gdb9_database.sh # becomes $0 inside the script above; unused, required by `bash -c`

hints:
  DockerRequirement:
    dockerPull: ndonyapour/molgan

inputs:
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
      position: 3

  output_NP_Score_path:
    label: Output ceout file (AMBER ceout)
    doc: |-
      Output ceout file (AMBER ceout)
      Type: string
      File type: output
      Accepted formats: gz
      Example file: https://github.com/bioexcel/biobb_amber/raw/master/biobb_amber/test/data/cphstats/sander.ceout.gz
    type: string?
    format:
    - edam:format_3987
    default: NP.gz
    inputBinding:
      position: 4

  output_SA_Score_path:
    label: Output ceout file (AMBER ceout)
    doc: |-
      Output ceout file (AMBER ceout)
      Type: string
      File type: output
      Accepted formats: gz
      Example file: https://github.com/bioexcel/biobb_amber/raw/master/biobb_amber/test/data/cphstats/sander.ceout.gz
    type: string?
    format:
    - edam:format_3987
    default: SA.gz
    inputBinding:
      position: 5


outputs:
  output_sdf_path:
    label: Path to the output file
    doc: |-
      Path to the output file
    type: File
    outputBinding:
      glob: $(inputs.output_sdf_path)
    format: edam:format_3814 # sdf

  output_NP_Score_path:
    label: Output ceout file 
    doc: |-
      Output ceout file 
    type: File
    outputBinding:
      glob: $(inputs.output_NP_Score_path)
    format: edam:format_3987 # gz

  output_SA_Score_path:
    label: Output ceout file 
    doc: |-
      Output ceout file 
    type: File
    outputBinding:
      glob: $(inputs.output_SA_Score_path)
    format: edam:format_3987 # gz

$namespaces:
  edam: https://edamontology.org/

$schemas:
- https://raw.githubusercontent.com/edamontology/edamontology/master/EDAM_dev.owl
