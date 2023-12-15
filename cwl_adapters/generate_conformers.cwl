#!/usr/bin/env cwl-runner
cwlVersion: v1.0

class: CommandLineTool

label: Download the PDBbind refined database

doc: |-
  Download the PDBbind refined database

baseCommand: ['python3', '/generate_conformers.py']

hints:
  DockerRequirement:
    dockerPull: ndonyapour/generate_conformers

requirements:
  InlineJavascriptRequirement: {}

inputs:
  input_excel_path:
    label: Path to the input xlsx file
    type: File
    format:
    - edam:format_3620
    inputBinding:
      prefix: --input_excel_path

  query:
    label: query str to search the dataset
    doc: |-
      query str to search the dataset
      Type: string
      File type: input
      Accepted formats: txt
    type: string
    format:
    - edam:format_2330
    inputBinding:
      prefix: --query

  output_txt_path:
    label: Path to the text dataset file
    doc: |-
      Path to the text dataset file
      Type: string
      File type: output
      Accepted formats: txt
    type: string
    format:
    - edam:format_2330
    inputBinding:
      prefix: --output_txt_path
    default: system.log

  output_sdf_path:
    label: Path to the input file
    doc: |-
      Path to the input file
      Type: string
      File type: input
      Accepted formats: sdf
    type: string
    format:
    - edam:format_3814 # sdf

  min_row:
    label: The row min index
    doc: |-
      The row min inex
      Type: int
    type: int?
    format:
    - edam:format_2330
    inputBinding:
      prefix: --min_row

  max_row:
    label: The row max index
    doc: |-
      The row max inex
      Type: int
    type: int?
    format:
    - edam:format_2330
    inputBinding:
      prefix: --max_row

  smiles_column:
    label: The name of the smiles column
    doc: |-
      The name of the smiles column
      Type: string
      File type: output
      Accepted formats: txt
    type: string
    format:
    - edam:format_2330
    inputBinding:
      prefix: --smiles_column

  binding_data_column:
    label: The name of the binding data column
    doc: |-
      The name of the binding data column
      Type: string
      File type: output
      Accepted formats: txt
    type: string
    format:
    - edam:format_2330
    inputBinding:
      prefix: --binding_data_column

  convert_Kd_dG:
    label: If this is set to true, dG will be calculated
    doc: If this is set to true, dG will be calculated  
    type: boolean
    format:
    - edam:format_2330
    inputBinding:
      prefix: --convert_Kd_dG
    default: False

  experimental_dGs:
    label: Experimental Free Energies of Binding
    doc: |-
      Experimental Free Energies of Binding
    type: string?
    format:
    - edam:format_2330

outputs:
  output_txt_path:
    label: Path to the txt file
    doc: |-
      Path to the txt file
    type: File
    outputBinding:
      glob: $(inputs.output_txt_path)
    format: edam:format_2330

  output_sdf_path:
    label: Path to the input file
    doc: |-
      Path to the input file
      Type: string
      File type: input
      Accepted formats: sdf
    type: File[]
    outputBinding:
      # NOTE: Do NOT just use glob: ./*.sdf !!! This will return an array sorted by filenames.
      # We want the order of output_sdf_paths to match the order of experimental_dGs, etc
      # Because we need to compare experimental ΔGs with predicted values.
      glob: $(inputs.output_txt_path)
      loadContents: true
      outputEval: |
        ${
          var lines = self[0].contents.split("\n");
          var sdfs = [];
          for (var idx = 0; idx < lines.length; idx++) {
            var words = lines[idx].split(" ");
            var sdffile = {"class": "File", "path": "ligand_" + idx + ".sdf"};
            sdfs.push(sdffile);
            }
            
          return sdfs;
        }
    format: edam:format_3814

  experimental_dGs:
    label: Experimental Free Energies of Binding
    doc: |-
      Experimental Free Energies of Binding
    type: float[]
    outputBinding:
      # NOTE: Do NOT just use $(inputs.output_txt_path) !!! This will return an array sorted by filenames.
      # We want the order of output_sdf_paths to match the order of experimental_dGs, etc
      # Because we need to compare experimental ΔGs with predicted values.
      glob: $(inputs.output_txt_path)
      loadContents: true
      outputEval: |
        ${
          var lines = self[0].contents.split("\n");
          var experimental_dGs = [];
          for (var i = 0; i < lines.length; i++) {
            var words = lines[i].split(" ");
            if (words.length > 2) {
              var experimental_dG = parseFloat(words[2]);
              experimental_dGs.push(experimental_dG);
            }
          }

          if (experimental_dGs.length == 0) {
            throw new Error("Error! Experimental dGs are empty!");
          } else {
            return experimental_dGs;
          }
        }

$namespaces:
  edam: https://edamontology.org/

$schemas:
- https://raw.githubusercontent.com/edamontology/edamontology/master/EDAM_dev.owl