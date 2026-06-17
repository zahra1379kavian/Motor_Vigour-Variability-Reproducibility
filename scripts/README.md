# Scripts

These are the files you run first. They have simple figure names and call the
longer code in `analysis/`.

Run scripts from the repository root. For example:

```bash
python scripts/figure_03_vigour_network_map.py
python scripts/supp_figure_09_medication_fc_consistency.py
```

Some entry points can use the packaged derived maps or processed result tables.
Others require subject-level beta, behavioural, BOLD, atlas, or intermediate
analysis inputs that are not included in this curated repository. Those scripts
will report missing inputs unless the external paths described in
`data/external/README.md` are available or replaced with command-line
arguments.

`supp_figure_04_standard_glm_feat.py` checks that the curated FEAT figure and
provenance files are present. Re-running that FSL analysis requires a working
FSL installation and the external subject-level BOLD inputs.
