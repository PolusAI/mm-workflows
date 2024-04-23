# Subworkflows e.g. "macros"

In the previous tutorial, we listed all of the workflow steps in a single file. Alternatively, we can extract some of the steps into another workflow.

<table>
<tr>
<td>
docs/tutorials/multistep3.wic

```yaml
steps:
- pdb:
    in:
      config: !ii
        pdb_code: 1aki
- fix_protein.wic:
- extract_model:
    in:
      config: !ii
        models: [1]
      output_structure_path: !ii protein_model_1.pdb
```

docs/tutorials/fix_protein.wic

```yaml
steps:
- fix_amides:
- fix_side_chain:
```

</td>
<td>
docs/tutorials/multistep3.wic.gv.png

![Multistep](multistep3.wic.gv.png)

</td>
</tr>
</table>

We have simply moved the fix_* steps into `fix_protein.wic` and called it from the main workflow. As you can see from the arrows in the graphical representation, the exact same edges have been inferred! The inference algorithm is guaranteed to work identically across subworkflow boundaries! You are completely free to abstract away minor details behind a subworkflow, and the main workflow graph will be identical.