#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

label: OnionNet (version1) for feature generation of docking poses

baseCommand: ["python", "/onionnet/generate_features.py"]

hints:
  DockerRequirement:
    dockerPull: cyangnyu/onionnet

requirements:
  InlineJavascriptRequirement: {}

inputs:
  complex_path_file:
    label: path file of protein-ligand complexes (structures in pdb format)
    type: File?
    format:
    - edam:format_1476
    inputBinding:
      prefix: -inp

  num_of_cpus:
    label: number of CPUs to use.
    type: int?
    format:
    - edam:format_2330
    inputBinding:
      prefix: -nt
    default: 1

  output_feature_file:
    label: the output file name containing the features.
    type: string?
    format:
    - edam:format_3752
    inputBinding:
      prefix: -out
    default: "output.csv"

outputs:
  output_feature_file:
    type: File
    format: edam:format_3752
    outputBinding:
      glob: $(inputs.output_feature_file)

$namespaces:
  edam: https://edamontology.org/
  cwltool: http://commonwl.org/cwltool#

$schemas:
- https://raw.githubusercontent.com/edamontology/edamontology/master/EDAM_dev.owl
