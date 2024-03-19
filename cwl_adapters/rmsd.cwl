cwlVersion: v1.2

class: CommandLineTool

baseCommand: ["python3", "/rmsd.py"]

hints:
  DockerRequirement:
    dockerPull: mrbrandonwalker/rmsd

inputs:

  predicted_poses:
    type: File[]
    inputBinding:
        itemSeparator: ","
        position: 1
        prefix: "--predicted_poses"

  reference_pose:
    type: File
    inputBinding:
        position: 2
        prefix: "--reference_pose"

  output_json_name:
    type: string
    default: "output.json"
    format:
    - edam:format_3816
    inputBinding:
      position: 3
      prefix: --output_json_name

  rmsd_output:
    type: string?

outputs:

  rmsd_output:
    type: File
    outputBinding:
      glob: $(inputs.output_json_name)

