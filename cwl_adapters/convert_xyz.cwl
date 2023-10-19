#!/usr/bin/env cwl-runner
cwlVersion: v1.2

class: CommandLineTool

label: This class is a wrapper of the Open Babel tool.

doc: |-
  Converts small molecules in 2D or 3D formats to the XYZ format.  

baseCommand: obabel
#Usage:
#obabel[-i<input-type>] <infilename> [-o<output-type>] -O<outfilename> [Options]
#...
#Options, other than -i -o -O -m, must come after the input files.
arguments: [$(inputs.input_path), "-o", "xyz", "-O", $(inputs.output_xyz_path)]

hints:
  DockerRequirement:
    dockerPull: quay.io/biocontainers/biobb_chemistry:4.0.0--pyhdfd78af_1

requirements:
  InlineJavascriptRequirement: {}

inputs:
  input_path:
    label: Path to the input file
    doc: |-
      Path to the input file
      Type: string
      File type: input
      Accepted formats: gro, mol, mol2, pdb, pdbqt, sdf
      Example file: https://github.com/bioexcel/biobb_chemistry/raw/master/biobb_chemistry/test/data/babel/babel.smi
    type: File
    format:
    - edam:format_2033
    - edam:format_3815
    - edam:format_3816
    - edam:format_1476
    - edam:format_1476
    - edam:format_3814

  output_xyz_path:
    label: Path to the output file
    doc: |-
      Path to the output file
      Type: string
      File type: output
      Accepted formats: xyz
    type: string
    format:
    - edam:format_3877 # xyz
    default: system.xyz

outputs:
  output_xyz_path:
    label: Path to the output file
    doc: |-
      Path to the output file
    type: File
    outputBinding:
      glob: $(inputs.output_xyz_path)
    format: edam:format_3877 # xyz

$namespaces:
  edam: https://edamontology.org/

$schemas:
- https://raw.githubusercontent.com/edamontology/edamontology/master/EDAM_dev.owl
