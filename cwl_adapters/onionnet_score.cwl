#!/usr/bin/env cwl-runner

cwlVersion: v1.0
class: CommandLineTool

label: OnionNet (version1) for rescoring of docking poses

baseCommand: ["conda", "run", "-n", "py37", "python", "/onionnet/predict.py"]

hints:
  DockerRequirement:
    dockerPull: ndonyapour/onionnet

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
    type: float[]
    outputBinding:
      glob: $(inputs.output_score_file)
      loadContents: true
      outputEval: |
        ${
          const lines = self[0].contents.split("\n");
          // Remove blank lines
          var non_blank_lines = lines.filter(function(line) {return line.trim() !== '';});
          var onionnet_dGs = [];   

          // Calculate the binding free energy from pKa, 
          // See https://github.com/PolusAI/mm-workflows/blob/1c2284e5bf517b90950c2f7af574caff3aa536d8/examples/scripts/generate_conformers.py#L31
          // Calculate pKa = -log(Kd), See https://en.wikipedia.org/wiki/Dissociation_constant,
          // Calculate the binding free energy from Kd so we can make the correlation plots. dG = RT*ln(Kd/c)
          // See https://en.wikipedia.org/wiki/Binding_constant
          const ideal_gas_constant = 8.31446261815324;  // J/(Mol*K)
          const kcal_per_joule = 4184;
          // NOTE: Unfortunately, the temperature at which experimental Kd binding data was taken
          // is often not recorded. Thus, we are forced to guess. The two standard guesses are
          // physiological body temperature (310K) or room temperature (298K).
          const temperature = 298;
          const RT = (ideal_gas_constant / kcal_per_joule) * temperature;
          // NOTE: For performance, simulations are often done in a very small unit cell, and
          // thus at a very high concentration. The size of the unit cell bounds the volume.
          // For shorter simulations where the ligand has not explored the entire box, it may
          // be less. See the Yank paper for a method of calculating the correct volumes.
          const standard_concentration = 1;  // Units of mol / L, but see comment above.  

          for (var i = 1; i < non_blank_lines.length; i++) {
            // The correct lines of output score file should be of the form
            // pKa_predicted 
            // /var/lib/cwl/stg19c300d1-f7fd-4a38-80d2-0f5615e3eb8f/complex_pdbs.pdb,7.441
            var pKa_string = non_blank_lines[i].split(",").filter(function(s) {return !isNaN(parseFloat(s))})[0];
            var Kd = (10 ** - parseFloat(pKa_string));
            var dG = RT * Math.log(Kd / standard_concentration);
            onionnet_dGs.push(dG);
          }
          return onionnet_dGs;
        }

$namespaces:
  edam: https://edamontology.org/
  cwltool: http://commonwl.org/cwltool#

$schemas:
- https://raw.githubusercontent.com/edamontology/edamontology/master/EDAM_dev.owl
