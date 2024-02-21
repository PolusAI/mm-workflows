#!/usr/bin/env cwl-runner
cwlVersion: v1.0

class: CommandLineTool

label: Download the PDBbind refined database

doc: |-
  Download the PDBbind refined database

baseCommand: python3

hints:
  DockerRequirement:
    dockerImageId: pdbbind_refined_v2020  # NOTE: no username
    dockerFile:
        $include: ../../../mm-workflows/examples/scripts/Dockerfile_pdbbind_refined

requirements:
  InlineJavascriptRequirement: {}

inputs:
  script:
    type: string
    inputBinding:
      position: 1
    default: /generate_pdbbind_complex.py

  index_file_name:
    label: The index file name
    type: string
    format:
    - edam:format_2330
    inputBinding:
      prefix: --index_file_name
      position: 2
    default: INDEX_refined_data.2020

  base_dir:
    label: The base_dir path
    type: string
    format:
    - edam:format_2330
    inputBinding:
      prefix: --base_dir
      position: 3
    default: /refined-set

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
      position: 4

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
      position: 5
    default: system.log

  output_pdb_paths:
    label: Path to the input file
    doc: |-
      Path to the input file
      Type: string
      File type: input
      Accepted formats: pdb
    type: string
    format:
    - edam:format_1476 # pdb

  output_sdf_paths:
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
      position: 6
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
      position: 7
      prefix: --max_row

  convert_Kd_dG:
    label: If this is set to true, dG will be calculated
    doc: If this is set to true, dG will be calculated
    type: boolean
    format:
    - edam:format_2330
    inputBinding:
      prefix: --convert_Kd_dG
      position: 8
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

  output_pdb_paths:
    label: Path to the input file
    doc: |-
      Path to the input file
      Type: string
      File type: input
      Accepted formats: pdb
    type: File[]
    outputBinding:
      # NOTE: Do NOT just use glob: ./*.pdb !!! This will return an array sorted by filenames.
      # We want the order of output_pdb_paths to match the order of experimental_dGs, etc
      # Becasue we need to compare experimental ΔGs with predicted values.
      glob: $(inputs.output_txt_path)
      loadContents: true
      outputEval: |
        ${
          var lines = self[0].contents.split("\n");
          var pdbs = [];
          for (var i = 0; i < lines.length; i++) {
            var words = lines[i].split(" ");
            var pdbid = words[0];
            var pdbfile = {"class": "File", "path": pdbid + "_protein.pdb"};
            pdbs.push(pdbfile);
            }

          return pdbs;
        }
    format: edam:format_1476

  output_sdf_paths:
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
          for (var i = 0; i < lines.length; i++) {
            var words = lines[i].split(" ");
            var pdbid = words[0];
            var sdffile = {"class": "File", "path": pdbid + "_ligand.sdf"};
            sdfs.push(sdffile);
            }
            
          return sdfs;
        }
    format: edam:format_3814

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
$namespaces:
  edam: https://edamontology.org/

$schemas:
- https://raw.githubusercontent.com/edamontology/edamontology/master/EDAM_dev.owl