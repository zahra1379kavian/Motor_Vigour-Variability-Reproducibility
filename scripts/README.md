# Scripts

These are short figure-facing entry points. Run them from the repository root:

```bash
python scripts/figure_03_vigour_network_map.py
python scripts/figure_04_ablation_anatomy.py
python scripts/figure_05_ablation_roi_summary.py
python scripts/supp_figure_09_medication_fc_consistency.py
```

Scripts that depend on subject-level data will report missing inputs if the
external data paths in `data/external/README.md` are not available.

Supplementary Figure 4 has no script yet because the requested FEAT output was
not present in the source repository.
