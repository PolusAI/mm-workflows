
#!/usr/bin/env cwl-runner
cwlVersion: v1.1  # See `loadContents: true` below!

class: CommandLineTool

label: Wrapper of the AutoDock Vina software.

doc: |-
  This class performs docking of the ligand to a set of grids describing the target protein via the AutoDock Vina software.

baseCommand: vina # NOTE: Only version >=1.2 supports --batch!
arguments:
# Need to parse box.pdb and pass in each number separately.
# REMARK BOX CENTER:     0.102     0.019    -0.004 SIZE:    30.195    31.940    27.005
# - "--dir" # Need to explicitly pass --dir . in --batch mode
# - "."
- "--center_x"
- $(inputs.input_box_path.contents.split("\n")[0].split(" ").filter(function(s) {return !isNaN(parseFloat(s))})[0])
- "--center_y"
- $(inputs.input_box_path.contents.split("\n")[0].split(" ").filter(function(s) {return !isNaN(parseFloat(s))})[1])
- "--center_z"
- $(inputs.input_box_path.contents.split("\n")[0].split(" ").filter(function(s) {return !isNaN(parseFloat(s))})[2])
- "--size_x"
- $(inputs.input_box_path.contents.split("\n")[0].split(" ").filter(function(s) {return !isNaN(parseFloat(s))})[3])
- "--size_y"
- $(inputs.input_box_path.contents.split("\n")[0].split(" ").filter(function(s) {return !isNaN(parseFloat(s))})[4])
- "--size_z"
- $(inputs.input_box_path.contents.split("\n")[0].split(" ").filter(function(s) {return !isNaN(parseFloat(s))})[5])
# NOTE: Cannot use a single javascript expression to create the entire arguments list because CWL treats it as a string:
# "the `arguments` field is not valid because value is a str"
#  ${
#    var words = inputs.input_box_path.contents.split("\n")[0].split(" ");
#    var nums = words.filter(function(s) {return !isNaN(parseFloat(s))});
#    var args = {};
#    args.push("--dir"); // Need to explicitly pass --dir . in --batch mode
#    args.push(".");
#    args.push("--center_x");
#    args.push(nums[0]);
#    args.push("--center_y");
#    args.push(nums[1]);
#    args.push("--center_z");
#    args.push(nums[2]);
#    args.push("--size_x");
#    args.push(nums[3]);
#    args.push("--size_y");
#    args.push(nums[4]);
#    args.push("--size_z");
#    args.push(nums[5]);
#    return args;
#  }

hints:
  DockerRequirement:
    dockerPull: jakefennick/autodock_vina

requirements:
  InlineJavascriptRequirement: {}

inputs:
  input_ligand_pdbqt_path:
    label: Path to the input PDBQT ligand
    doc: |-
      Path to the input PDBQT ligand
      Type: string
      File type: input
      Accepted formats: pdbqt
      Example file: https://github.com/bioexcel/biobb_vs/raw/master/biobb_vs/test/data/vina/vina_ligand.pdbqt
    type: File
    format:
    - edam:format_1476
    inputBinding:
      prefix: --ligand

  input_receptor_pdbqt_path:
    label: Path to the input PDBQT receptor
    doc: |-
      Path to the input PDBQT receptor
      Type: string
      File type: input
      Accepted formats: pdbqt
      Example file: https://github.com/bioexcel/biobb_vs/raw/master/biobb_vs/test/data/vina/vina_receptor.pdbqt
    type: File
    format:
    - edam:format_1476
    inputBinding:
      #position: 2
      prefix: --receptor

  input_box_path:
    label: Path to the PDB containing the residues belonging to the binding site
    doc: |-
      Path to the PDB containing the residues belonging to the binding site
      Type: string
      File type: input
      Accepted formats: pdb
      Example file: https://github.com/bioexcel/biobb_vs/raw/master/biobb_vs/test/data/vina/vina_box.pdb
    type: File
    format:
    - edam:format_1476
#    inputBinding:
#      position: 3
#      prefix: --input_box_path
    loadContents: true  # requires cwlVersion: v1.1
    # See https://www.commonwl.org/v1.1/CommandLineTool.html#Changelog
    # Curiously, cwlVersion: v1.0 allows loadContents for outputs, but not inputs.

  output_pdbqt_path:
    label: Path to the output PDBQT file
    doc: |-
      Path to the output PDBQT file
      Type: string
      File type: output
      Accepted formats: pdbqt
      Example file: https://github.com/bioexcel/biobb_vs/raw/master/biobb_vs/test/reference/vina/ref_output_vina.pdbqt
    type: string
    format:
    - edam:format_1476
    inputBinding:
      prefix: --out
    default: system.pdbqt

  output_log_path:
    label: Path to the log file
    doc: |-
      Path to the log file
      Type: string
      File type: output
      Accepted formats: log
      Example file: https://github.com/bioexcel/biobb_vs/raw/master/biobb_vs/test/reference/vina/ref_output_vina.log
    type: string
    format:
    - edam:format_2330
    default: system.log

  cpu:
    label: The number of CPUs to use
    doc: |-
      The number of CPUs to use
      Type: int
    type: int?
    format:
    - edam:format_2330
    inputBinding:
      prefix: --cpu
    default: 1

  exhaustiveness:
    label: exhaustiveness of the global search (roughly proportional to time).
    doc: |-
      exhaustiveness of the global search (roughly proportional to time).
      Type: int
    type: int?
    format:
    - edam:format_2330
    inputBinding:
      prefix: --exhaustiveness
    default: 8

  num_modes:
    label: maximum number of binding modes to generate
    doc: |-
      Tmaximum number of binding modes to generate
      Type: int
    type: int?
    format:
    - edam:format_2330
    inputBinding:
      prefix: --num_modes
    default: 9 

  min_rmsd:
    label: The minimum RMSD between output poses     
    doc: |-
      The minimum RMSD between output poses
      Type: int
    type: int?
    format:
    - edam:format_2330
    inputBinding:
      prefix: --min_rmsd
    default: 1

  energy_range:
    label: maximum energy difference between the best binding mode and the worst one displayed (kcal/mol).
    doc: |-
      maximum energy difference between the best binding mode and the worst one displayed (kcal/mol).
      Type: int
    type: int?
    format:
    - edam:format_2330
    inputBinding:
      prefix: --energy_range
    default: 3

outputs:
  output_pdbqt_path:
    label: Path to the output PDBQT file
    doc: |-
      Path to the output PDBQT file
    type: File
    outputBinding:
      glob: $(inputs.output_pdbqt_path)
    format: edam:format_1476

  output_log_path:
    label: Path to the log file
    doc: |-
      Path to the log file
    type: File
    outputBinding:
      glob: $(inputs.output_log_path)
    format: edam:format_2330

stdout: $(inputs.output_log_path)

$namespaces:
  edam: https://edamontology.org/

$schemas:
- https://raw.githubusercontent.com/edamontology/edamontology/master/EDAM_dev.owl
