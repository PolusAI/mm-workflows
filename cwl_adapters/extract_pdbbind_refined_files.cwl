#!/usr/bin/env cwl-runner
cwlVersion: v1.0

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
  InitialWorkDirRequirement:
    listing: |
      ${
        var lst = [];
        // console.log("inputs.data_dir");
        // console.log(inputs.data_dir)
        console.log("inputs.data_dir.listing.length.toString()");
        console.log(inputs.data_dir.listing.length.toString());
        for (var i = 0; i < inputs.data_dir.listing.length; i++) {
          // console.log(inputs.data_dir.listing[i].basename);

          if (inputs.pdb_ids.includes(inputs.data_dir.listing[i].basename)) {
            // console.log(inputs.data_dir.listing[i]);
            var dict = {
              "entry": inputs.data_dir.listing[i],
              "writable": true // Important!
            };
            lst.push(dict);
          }
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