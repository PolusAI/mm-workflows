#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

label: OnionNet (version1) for rescoring of docking poses

baseCommand: ["python", "/onionnet/predict.py"]

hints:
  DockerRequirement:
    dockerPull: cyangnyu/onionnet

requirements:
  InlineJavascriptRequirement: {}

inputs:
  input_feature_file:
    label: feature csv file for protein-ligand complexes
    type: File?
    format:
    - edam:format_3752
    inputBinding:
      prefix: -fn

  scaler:
    label: the standard scaler file.
    type: string?
    format:
    - edam:format_2330
    inputBinding:
      prefix: -scaler
    default: "/onionnet/models/StandardScaler.model"

  weights:
    label: the trained DNN model file.
    type: string?
    format:
    - edam:format_2330
    inputBinding:
      prefix: -weights
    default: "/onionnet/models/CNN_final_model_weights.h5"

  output_score_file:
    label: the predicted pKa values file
    type: string?
    format:
    - edam:format_3752
    inputBinding:
      prefix: -out
    default: "predicted_pKa.csv"

  onionnet_score:
    type: string?

outputs:
  output_score_file:
    type: File
    outputBinding:
      glob: $(inputs.output_score_file)
    format: edam:format_3752

  onionnet_score:
    label: Estimated Free Energy of Binding (onionnet score)
    doc: |-
      Estimated Free Energy of Binding
    type: float
    outputBinding:
      glob: $(inputs.output_score_file)
      loadContents: true
      outputEval: |
        ${
          const lines = self[0].contents.split("\n");
          // The correct line should be of the form
          // ,pKa_predicted
          // /var/lib/cwl/stg19c300d1-f7fd-4a38-80d2-0f5615e3eb8f/complex_pdbs.pdb,7.441
          const bfe_line = lines[1];
          // refactor can be used to convert pKa to binding free enegy, based on deltaG = -RT*lnK 
          const refactor = -0.73349;
          const docking_score_string = bfe_line.split(",").filter(function(s) {return !isNaN(parseFloat(s))})[0];
          const onionnet_score = parseFloat(docking_score_string)/refactor;
          return onionnet_score
        }

$namespaces:
  edam: https://edamontology.org/
  cwltool: http://commonwl.org/cwltool#

$schemas:
- https://raw.githubusercontent.com/edamontology/edamontology/master/EDAM_dev.owl
