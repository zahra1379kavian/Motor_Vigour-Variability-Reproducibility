# Figure 8 GVS Effects

This directory contains the companion result files for Figure 8.

- `rt_lme/`: trial table, mixed-model summaries, pairwise GVS-vs-sham tests,
  and the original Figure 8A export.
- `vigour_network_feature_delta/`: GVS feature summaries for the vigour-network
  BOLD projection and the original Figure 8B export.
- `task_activation_feature_delta/`: GVS feature summaries for the task-activation
  BOLD signal and the original Figure 8C export.

Run `python scripts/figure_08_gvs_effects.py` from the repository root to
refresh the curated panels in `figures/main/`. The script replots Figure 8A from
the packaged trial table and preserves the exact curated Figure 8B/8C exports by
copying them from this result bundle. Direct module users can pass
`--replot-feature-panels` to regenerate Figure 8B/8C from the packaged CSV
tables.
