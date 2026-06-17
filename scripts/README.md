# Scripts

These are the files you run first. They have simple figure names and call the
longer code in `analysis/`.

Run them from the repository root:

```bash
python scripts/figure_03_vigour_network_map.py
python scripts/figure_04_ablation_anatomy.py
python scripts/figure_05_ablation_roi_summary.py
python scripts/supp_figure_09_medication_fc_consistency.py
```

Scripts that depend on subject-level data will report missing inputs if the
external data paths in `data/external/README.md` are not available.

`supp_figure_04_standard_glm_feat.py` checks that the curated FEAT figure and
provenance files are present. Re-running that FSL analysis requires a working
FSL installation and the external subject-level BOLD inputs.
