#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

label: Smina docking tool to perform protein-ligand docking

baseCommand: ["smina"]

hints:
  DockerRequirement:
    dockerPull: cyangnyu/smina

requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - $(inputs.ligand_file)

inputs:

  receptor_file:
    label: receptor structure.
    doc: |-
      Path to the input receptor file
      Type: string?
      File type: input
      Accepted formats: pdb, pdbqt, sdf, mol, mol2
    type: File?
    format:
    - edam:format_1476
    - edam:format_3814
    - edam:format_3815
    - edam:format_3816
    inputBinding:
      prefix: -r

  ligand_file:
    label: ligand structure.
    doc: |-
      Path to the input ligand file
      Type: string?
      File type: input
      Accepted formats: pdb, pdbqt, sdf, mol, mol2
    type: File?
    format:
    - edam:format_1476
    - edam:format_3814
    - edam:format_3815
    - edam:format_3816
    inputBinding:
      prefix: -l

  ligand_box:
    label: 3D structure to define ligand box center and size
    doc: |-
      Path to the input docking box file
      Type: string?
      File type: input
      Accepted formats: pdb, pdbqt, sdf, mol, mol2
    type: File?
    format:
    - edam:format_1476
    - edam:format_3814
    - edam:format_3815
    - edam:format_3816
    inputBinding:
      prefix: --autobox_ligand

  local_only:
    label: try local minimization only rather than docking
    type: boolean?
    inputBinding:
      prefix: --local_only

  score_only:
    label: Do not do any conformational search; simply rescore.
    type: boolean?
    inputBinding:
      prefix: --score_only

  scoring:
    label: scoring function option, default is vina, options can be (vina, vinardo, or a customized scoring function) 
    type: string?
    inputBinding:
      prefix: --scoring
    default: "vina"

  output_dock_file:
    label: output docking poses.
    type: string?
    format: edam:format_1476
    inputBinding:
      prefix: -o
    default: "docked.pdb"

# Is this still necessary?
  output_path:
    type: string
    
outputs:
  
  output_dock_file:
    type: File
    outputBinding:
      glob: $(inputs.output_dock_file)
    format: edam:format_1476

  output_path:
    type: File
    outputBinding:
      glob: $(inputs.output_path)
    

stdout: $(inputs.output_path)

$namespaces:
  edam: https://edamontology.org/
  cwltool: http://commonwl.org/cwltool#

$schemas:
- https://raw.githubusercontent.com/edamontology/edamontology/master/EDAM_dev.owl
