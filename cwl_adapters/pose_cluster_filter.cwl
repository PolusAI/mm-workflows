#!/usr/bin/env cwl-runner
cwlVersion: v1.0

class: CommandLineTool

label: Cluster poses in protein and take max confident pose for each cluster

doc: |-
  Cluster poses in protein and take max confident pose for each cluster

baseCommand: ["python", "/pose_cluster_filter.py"]

hints:
  DockerRequirement:
    dockerPull: mrbrandonwalker/pose_cluster_filter

requirements:
  InlineJavascriptRequirement: {}

inputs:

  centroid_cutoff:
    type: float
    label: if centroid of all poses are within cutoff then keep only most confident pose
    default: 5
    inputBinding:
      prefix: --centroid_cutoff

  predicted_poses:
    type: File[]
    inputBinding:
      prefix: --predicted_poses
    format:
    - edam:format_3814

outputs:

  filtered_poses:
    type: File[]
    label: filtered poses
    outputBinding:
      glob: filtered_poses.txt # This determines what binds to self[0]
      loadContents: true
      outputEval: |
        ${
          // file contents look like
          // file_index cluster_index
          const lines = self[0].contents.split("\n").filter(line => line.trim() !== '');
          const lst = [];
          for (var i = 0; i < lines.length; i++) {
            var splitLine = lines[i].split(" ");
            // now find the File from inputs.predicted_poses
            var mol_idx = parseInt(splitLine[0]);
            var file = inputs.predicted_poses[mol_idx];
            lst.push(file);
          }
          return lst;
        }
    format: edam:format_3814

$namespaces:
  edam: https://edamontology.org/

$schemas:
- https://raw.githubusercontent.com/edamontology/edamontology/master/EDAM_dev.owl