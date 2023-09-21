#!/usr/bin/env cwl-runner
cwlVersion: v1.0

class: CommandLineTool

label: Downloads the NCI lgand database(s)

doc: |-
  Downloads the NCI lgand database(s)

# cwltool uses the command, arguments, resolved paths of inputs of File type, and
# a few requirement field (e.g. DockerRequirement, InitialWorkDirRequirement) to
# construct a keydict and calculate its md5 hashcode as the name of the cache
# folder of a workflow step. In order to reuse the cached output when rerunning,
# we should avoid using nondeterministic values, like $(runtime.outdir), and
# favor static values, like "." for output directory.
# See https://github.com/common-workflow-language/cwltool/blob/20f01e04328537714e57d136e242d3e7a9d44266/cwltool/command_line_tool.py#L850C1-L906C1.

baseCommand: cp
arguments: [$(inputs.database), $(".")]

requirements:
  InlineJavascriptRequirement: {}
  DockerRequirement:
    dockerPull: jakefennick/data

inputs:
  database:
    type: string
    default: /NCIOpen.sdf # NOTE: Initial / required
    # /NCIOpen.sdf is mounted inside the docker image, . is outside of the image.

outputs:
  output_sdf_path:
    label: Path to the output sdf file
    doc: |-
      Path to the output sdf file
    type: File
    format: edam:format_3814 # sdf
    outputBinding:
      glob: $(inputs.database.slice(1)) # Remove initial /

stdout: stdout

$namespaces:
  edam: https://edamontology.org/

$schemas:
- https://raw.githubusercontent.com/edamontology/edamontology/master/EDAM_dev.owl