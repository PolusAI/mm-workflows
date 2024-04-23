# Advanced Features

## Static dispatch

Here is an example that shows how to swap out constant pressure implementations.

```yaml
wic:
  default_implementation: gromacs
  implementations:
    gromacs:
      steps:
        - npt_gromacs.wic:
    amber:
      steps:
        - npt_amber.wic:
  graphviz:
    label: Constant Pressure
```

Then you just need to choose a specific implementation at the call site:

```yaml
steps:
  - nvt.wic:
  - npt.wic:

wic:
  graphviz:
    label: Equilibration
  steps:
    (2, npt.wic):
      wic:
        implementation: amber
```
This will override the default implementation of `gromacs` and use `amber`. This really just means that `npt_amber.wic` is called instead of `npt_gromacs.wic` (If `--insert_steps_automatically` is enabled, the compiler will attempt to automatically insert the necessary file format conversions as determined below.)

## Subinterpreters

A portion of [`examples/gromacs/nvt.wic`](https://github.com/PolusAI/mm-workflows/blob/main/examples/gromacs/nvt.wic) in `mm-workflows` is shown below. You can see that the `in:` tag of gmx_energy is identical to the `config:` tag of cwl_subinterpreter. This currently needs to be manually copy & pasted (and indented), but it should be possible to automatically do this in the future.

```yaml
...
  - mdrun:
      out:
      - output_edr_path: !& nvt.edr # Explicit edge reference / anchor
        # (This edge can be inferred, but made explicit for demonstration purposes.)
  - gmx_energy:
      in:
        input_energy_path: !* nvt.edr # Explicit edge dereference / alias
        config: !ii
          terms: [Temperature]
        output_xvg_path: temperature.xvg
# NOTE: explicit edges are not supported with cwl_subinterpreter, and all filenames
# must be globally unique!
  - cwl_subinterpreter:
      in:
        #cachedir_path: /absolute/path/to/cachedir/ (automatically filled in by wic)
        file_pattern: '*nvt.edr'  # Any strings that start with & or * need to be escaped in quotes
        cwl_tool: gmx_energy # This can also be an arbitrary subworkflow!
        max_times: '5'
        config: !ii
          in:
            input_energy_path: '*nvt.edr' # This * is automatically removed.
            config: !ii
              terms: [Temperature]
            output_xvg_path: temperature.xvg
...
```

Note that although gmx_energy appears before cwl_subinterpreter in the YAML file, gmx_energy is independent of cwl_subinterpreter in the DAG and thus not considered to be a previous step. We include gmx_energy simply to guarantee that the analysis gets run one more time in the main workflow, when all the files are known to be in their final state.

### Known Issues

Since the two runtimes are not linked, there is not currently a reliable way to determine if the previous steps have finished. Thus, to guarantee termination of the second runtime, we simply execute `cwl_tool` upto `max_times`. We also waive any guarantees about the files, so the subworkflow in the second runtime may of course fail for any number of reasons. Thus, we do not propagate speculative failures up to the main workflow.

The runtime system intentionally hides the working sub-directories of each step. Thus, we are forced to use a file watcher recursively starting from `cachedir_path`. This is why all filenames used with cwl_subinterpreter must be globally unique. (Actually, for technical reasons we cannot use a file watching library; we simply use a good old fashioned polling loop.)

## Real-time plots

It is assumed that the real-time analysis takes care of the complex log file parsing, etc and produces simple tabular data files (i.e. csv files separated by whitespace instead of a comma). We need to use the same file watching / polling trick as above to locate these tabular data files. The first argument to the following command is the directory in which to look for the files. (By default it is `cachedir` because that is the default value of the  `--cachedir` wic command line argument.) You can also optionally supply the file patterns, which by default are `*.xvg` and `*.dat`.

```
timeseriesplots cachedir <pat1> <pat2> <...>
```

## YAML Metadata Annotations

### Overloading / Parameter Passing

This example shows how we can recursively pass in parameters / recursively overload metadata.

Suppose we want to do a very careful minimization, first in vacuum and then in solvent (i.e. [`examples/gromacs/setup_vac_min.wic`](https://github.com/PolusAI/mm-workflows/blob/main/examples/gromacs/setup_vac_min.wic) in `mm-workflows`). We would like to re-use the abstract minimization protocol from `min.wic`. However, our stability analysis requires an explicit edge definition from the final minimized coordinates (i.e. in solvent). If we try to simply add `- output_tpr_path: !& min.tpr` directly to `min.wic`, there will be duplicate definitions! This is not allowed (it will generate an exception).

The solution is to pass in this parameter to only the second instance of `min.wic`.

A portion of [`examples/gromacs/basic.wic`](https://github.com/PolusAI/mm-workflows/blob/main/examples/gromacs/basic.wic) is shown below.

```yaml
...
# Put everything under one top-level wic: tag to facilitate easy merging and removal.
wic:
  graphviz:
    label: Molecular Dynamics
  steps:
    (1, min.wic):
      wic:
        steps:
          (2, cg.wic):
            wic:
              steps:
                (1, grompp):
                  out:
                  - output_tpr_path: !& min.tpr
...
```