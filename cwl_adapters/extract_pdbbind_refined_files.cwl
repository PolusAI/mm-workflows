#!/usr/bin/env cwl-runner
cwlVersion: v1.2

class: CommandLineTool

label: Download the PDBbind refined database

doc: |-
  Download the PDBbind refined database

baseCommand: echo

hints:
  cwltool:LoadListingRequirement:
    loadListing:
      no_listing:

requirements:
  InlineJavascriptRequirement: {}
  InplaceUpdateRequirement:
    inplaceUpdate: true
  InitialWorkDirRequirement:
    listing: |
      ${
        var lst = [];
        // console.log("inputs.data_dir");
        // console.log(inputs.data_dir)
        console.log("inputs.pdb_ids.length.toString()");
        console.log(inputs.pdb_ids.length.toString());
        for (var i = 0; i < inputs.pdb_ids.length; i++) {
          // console.log(inputs.pdb_ids[i]);
          var dict = {
            "class": "File",
            "location": inputs.data_dir.location.concat("/", inputs.pdb_ids[i]),
          };
          lst.push(dict);
        }
        console.log("lst.length.toString()");
        console.log(lst.length.toString());
        return lst;
      }

inputs:

  data_dir:
    type: Directory
    default:
      class: Directory
      location: ../../../mm-workflows/examples/scripts/refined-set

  pdb_ids:
    label: PDB ID's to extract files for
    type:
      type: array
      items: string


outputs:

  output_pdb_paths:
    label: Path to the PDB file
    format: edam:format_1476
    type: File[]
    outputBinding:
      glob: "*.pdb"

  output_sdf_paths:
    label: Path to the SDF file
    type: File[]
    outputBinding:
      glob: "*.sdf"
    format: edam:format_3814



$namespaces:
  edam: https://edamontology.org/
  cwltool: "http://commonwl.org/cwltool#"

$schemas:
- https://raw.githubusercontent.com/edamontology/edamontology/master/EDAM_dev.owl