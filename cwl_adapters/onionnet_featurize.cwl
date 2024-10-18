#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

label: Generate the features for the protein-ligand complexes using OnionNet V1

baseCommand: ["conda", "run", "-n", "py37", "python", "/onionnet/generate_features.py"]
arguments: ["-inp", "input.dat"]

hints:
  DockerRequirement:
    dockerPull: ndonyapour/onionnet

requirements:
  - class: InlineJavascriptRequirement
  # This enables staging input files and dynamically generating a file
  # containing the file paths on the fly.
  - class: InitialWorkDirRequirement 
    listing: |
      ${
        var dat_file_contents = "";
        for (var i = 0; i < inputs.pdb_paths.length; i++) {
          dat_file_contents += inputs.pdb_paths[i].path + "\n";
        }
        // Note: Uses https://www.commonwl.org/v1.0/CommandLineTool.html#Dirent
        // and https://www.commonwl.org/user_guide/topics/creating-files-at-runtime.html
        // "If the value is a string literal or an expression which evaluates to a string, 
        // a new file must be created with the string as the file contents."
        return ([{"entryname": "input.dat", "entry": dat_file_contents}]);
      }

inputs:
  pdb_paths:
    label: The path of input pdb files
    type: File[]
    format: 
    - edam:format_1476 # PDB

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
