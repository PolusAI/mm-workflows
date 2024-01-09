#!/usr/bin/env cwl-runner

# https://github.com/PolusAI/workflow-inference-compiler
steps:
  unzip_gro_amb__step__1__unzip_top:
    run: ../unzip_top.cwl
    in:
      input_top_zip_path: unzip_gro_amb__step__1__unzip_top__unzip_gro_amb__step__2__bss_gro_amb__output_top_path
    out:
    - output_top_path
  unzip_gro_amb__step__2__bss_gro_amb:
    run: bss_gro_amb.cwl
    in:
      input_top_path: unzip_gro_amb__step__1__unzip_top__output_top_path
      input_crd_path: unzip_gro_amb__step__2__bss_gro_amb__input_crd_path
    out:
    - output_top_path
    - output_crd_path
cwlVersion: v1.0
class: Workflow
$namespaces:
  edam: https://edamontology.org/
$schemas:
- https://raw.githubusercontent.com/edamontology/edamontology/master/EDAM_dev.owl
inputs:
  unzip_gro_amb__step__1__unzip_top__input_top_zip_path:
    type: File
    format:
    - edam:format_3987
  unzip_gro_amb__step__2__bss_gro_amb__input_crd_path:
    type: File
    format:
    - edam:format_2033
outputs:
  unzip_gro_amb__step__1__unzip_top__output_top_path:
    type: File
    format: edam:format_3880
    outputSource: unzip_gro_amb__step__1__unzip_top/output_top_path
  unzip_gro_amb__step__2__bss_gro_amb__output_top_path:
    type: File
    format: edam:format_3881
    outputSource: unzip_gro_amb__step__2__bss_gro_amb/output_top_path
  unzip_gro_amb__step__2__bss_gro_amb__output_crd_path:
    type: File
    format: edam:format_3878
    outputSource: unzip_gro_amb__step__2__bss_gro_amb/output_crd_path
