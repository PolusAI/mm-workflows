#!/usr/bin/env cwl-runner
cwlVersion: v1.0

class: CommandLineTool

label: Determine PDBbind refined database dG and filenames to extract

doc: |-
  Determine PDBbind refined database dG and filenames to extract

baseCommand: python3

hints:
  DockerRequirement:
    dockerPull: mrbrandonwalker/pdbbind_refined_v2020


requirements:
  InlineJavascriptRequirement: {}

inputs:

  script:
    type: string
    inputBinding:
      position: 1
    default: /generate_pdbbind_complex.py

  index_file_path:
    label: The index file path
    type: File
    inputBinding:
      prefix: --index_file_path
      position: 2
    default:
      class: File
      location: ../../../mm-workflows/examples/scripts/refined-set/index/INDEX_refined_data.2020

  query:
    label: query str to search the dataset, Pandas query doesn't support slash(/) in column names please use Kd_Ki instead of Kd/Ki
    doc: |-
      query str to search the dataset. Pandas query doesn't support slash(/) in column names please use Kd_Ki instead of Kd/Ki
      Type: string
      File type: input
      Accepted formats: txt
    type: string
    format:
    - edam:format_2330
    inputBinding:
      prefix: --query
      position: 3

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
      position: 4
    default: system.log

  min_row:
    label: The row min index
    doc: |-
      The row min inex
      Type: int
    type: int?
    format:
    - edam:format_2330
    inputBinding:
      position: 5
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
      position: 6
      prefix: --max_row

  convert_Kd_dG:
    label: If this is set to true, dG will be calculated
    doc: If this is set to true, dG will be calculated
    type: boolean
    format:
    - edam:format_2330
    inputBinding:
      prefix: --convert_Kd_dG
      position: 7
    default: False

  experimental_dGs:
    label: Experimental Free Energies of Binding
    doc: |-
      Experimental Free Energies of Binding
    type: string?
    format:
    - edam:format_2330

  pdb_ids:
    label: The PDBID of proteins
    doc: |-
      The PDBID of proteins
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

  experimental_dGs:
    label: Experimental Free Energies of Binding
    doc: |-
      Experimental Free Energies of Binding
    type: ["null", {"type": "array", "items": "float"}]
    outputBinding:
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
            return null;
          } else {
            return experimental_dGs;
          }
        }

  pdb_ids:
    label: The PDBID of proteins
    doc: |-
      The PDBID of proteins
    type:
      type: array
      items: string
    outputBinding:
      glob: $(inputs.output_txt_path)
      loadContents: true
      outputEval: |
        ${
          var lines = self[0].contents.split("\n");
          var pdbids = [];
          for (var i = 0; i < lines.length; i++) {
            var words = lines[i].split(" ");
            pdbids.push(words[0]);
            }

          if (pdbids.length == 0) {
            throw new Error("Error! pdbids are empty!");
          } else {
            return pdbids;
          }
        }

  stdout:
    type: File
    outputBinding:
      glob: stdout

stdout: stdout

$namespaces:
  edam: https://edamontology.org/

$schemas:
- https://raw.githubusercontent.com/edamontology/edamontology/master/EDAM_dev.owl