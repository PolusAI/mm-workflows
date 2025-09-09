#!/usr/bin/env cwl-runner
cwlVersion: v1.0

class: CommandLineTool

label: This class is a wrapper of the Open Babel tool.

doc: |-
  Small molecule format conversion for structures or trajectories. Open Babel is a chemical toolbox designed to speak the many languages of chemical data. It's an open, collaborative project allowing anyone to search, convert, analyze, or store data from molecular modeling, chemistry, solid-state materials, biochemistry, or related areas. Visit the official page.

baseCommand: obabel
#Usage:
#obabel[-i<input-type>] <infilename> [-o<output-type>] -O<outfilename> [Options]
#...
#Options, other than -i -o -O -m, must come after the input files.
arguments: [$(inputs.input_path), "-h", "-o", "pdb", "-O", $(inputs.output_pdb_path)]
# NOTE: These arguments must be given individually; they cannot be concatenated together.
# (e.g. -xrhn) Otherwise, all but the first argument will be silently ignored!

# pdb format options
# Read Options e.g. -as
#   s  Output single bonds only
#   b  Disable bonding entirely
#   c  Ignore CONECT records

# Write Options, e.g. -xo
#   n  Do not write duplicate CONECT records to indicate bond order
#   o  Write origin in space group label (CRYST1 section)

hints:
  DockerRequirement:
    dockerPull: quay.io/biocontainers/biobb_chemistry:4.0.0--pyhdfd78af_1

requirements:
  InlineJavascriptRequirement: {}

inputs:
  first_molecule:
    label: Index of the first molecule (1-based)
    doc: |-
      Input Index of the first molecule (1-based)
      Type: string?
    type: string?
    format:
    - edam:format_2330 # 'Textual format'
    inputBinding:
      prefix: -f

  last_molecule:
    label: Index of the last molecule (1-based)
    doc: |-
      Input Index of the last molecule (1-based)
      Type: string?
    type: string?
    format:
    - edam:format_2330 # 'Textual format'
    inputBinding:
      prefix: -l

  input_path:
    label: Path to the input file
    doc: |-
      Path to the input file
      Type: string
      File type: input
      Accepted formats: dat, ent, fa, fasta, gro, inp, log, mcif, mdl, mmcif, mol, mol2, pdb, pdbqt, png, sdf, smi, smiles, txt, xml, xtc
      Example file: https://github.com/bioexcel/biobb_chemistry/raw/master/biobb_chemistry/test/data/babel/babel.smi
    type: File
    format:
    - edam:format_1637
    - edam:format_1476
    - edam:format_1929
    - edam:format_1929
    - edam:format_2033
    - edam:format_3878
    - edam:format_2030
    - edam:format_1477
    - edam:format_3815
    - edam:format_1477
    - edam:format_3815
    - edam:format_3816
    - edam:format_1476
    - edam:format_1476
    - edam:format_3603
    - edam:format_3814
    - edam:format_1196
    - edam:format_1196
    - edam:format_2033
    - edam:format_2332
    - edam:format_3875

  output_pdb_path:
    label: Path to the output file
    doc: |-
      Path to the output file
      Type: string
      File type: output
      Accepted formats: pdb
    type: string
    format:
    - edam:format_1476 # pdb
    default: system.pdb

  arg1:
    label: Additional arguments
    doc: |-
      Additional arguments
      Type: string?
    type: string?
    format:
    - edam:format_2330 # 'Textual format'
    inputBinding:
      position: 1

outputs:
  output_pdb_path:
    label: Path to the output file
    doc: |-
      Path to the output file
    type: File
    outputBinding:
      glob: $(inputs.output_pdb_path)
    format: edam:format_1476 # pdb

$namespaces:
  edam: https://edamontology.org/

$schemas:
- https://raw.githubusercontent.com/edamontology/edamontology/master/EDAM_dev.owl
