# Multi-step workflows

Our first multi step workflow will consist of downloading a protein from an online database. Unfortunately, experiments are typically unable to resolve all of the atoms and/or residues, so it is necessary to 'fix' the initial data.

On the right is a visual representation of the workflow as a computational graph. The nodes are the steps, and the edges are the input and output files. This graph representation is automatically generated every time you compile a workflow (with `--graphviz`). It is very useful for visually debugging issues with workflows, and so it is a very good idea to ***`always look at the graph representation`*** before running a workflow.

<table>
<tr>
<td>
docs/tutorials/multistep1.wic

```yaml
steps:
- pdb:
    in:
      config: !ii
        pdb_code: 1aki
      output_pdb_path: !ii protein.pdb
    out:
    - output_pdb_path: !& protein_edge
- fix_amides:
    in:
      input_pdb_path: !* protein_edge
      output_pdb_path: !ii protein_fix_amides.pdb
    out:
    - output_pdb_path: !& protein_fix_amides_edge
- fix_side_chain:
    in:
      input_pdb_path: !* protein_fix_amides_edge
      output_pdb_path: !ii protein_fix_side_chain.pdb
    out:
    - output_pdb_path: !& protein_fix_side_chain_edge
- extract_model:
    in:
      config: !ii
        models: [1]
      input_structure_path: !* protein_fix_side_chain_edge
      output_structure_path: !ii protein_model_1.pdb
    out:
    - output_structure_path: !& protein_model_1_edge
```

</td>
<td>
docs/tutorials/multistep1.wic.gv.png

![Multistep](multistep1.wic.gv.png)

</td>
</tr>
</table>

## Explicit Edges

The first thing you might notice is the `!&` and `!*` notation. This is the syntax for explicitly creating an edge between an output and a later input. Note that `!&` must come first and can only be defined once, but then you can use `!*` multiple times in any later step.

## Edge Inference

Creating explicit edges can be a bit tedious and verbose, but in many cases the correct edges can be determined automatically. In this case, all of the steps take the previous pdb file as input and produce an output pdb file, so this is rather trivial. However, for more complex workflows edge inference can drastically simplify the wic files. Note that explicit edges are drawn in blue, and inferred edges are drawn in black/white.

For technical reasons edge inference is not perfect, so ***`users should always check that edge inference actually produces the intended graph`***.

<table>
<tr>
<td>
docs/tutorials/multistep2.wic

```yaml
steps:
- pdb:
    in:
      config: !ii
        pdb_code: 1aki
- fix_amides:
- fix_side_chain:
- extract_model:
    in:
      config: !ii
        models: [1]
      output_structure_path: !ii protein_model_1.pdb
```

</td>
<td>
docs/tutorials/multistep2.wic.gv.png

![Multistep](multistep2.wic.gv.png)

</td>
</tr>
</table>

Next, we will see how we can refactor some steps into a separate workflow.