# Intra-ROI vs Between-ROI FC Run-Baseline Comparison

This companion analysis uses the same ROI masks, voxel weights, and
intra-vs-between FC definitions as `med_effects.py`.

For each subject/session, FC was recomputed separately from beta run 1 and beta
run 2.

The plotted run-to-run baseline is a within-session variability magnitude:
`0.5 * (abs(OFF_run2 - OFF_run1) + abs(ON_run2 - ON_run1))`.

The plotted medication effect is recomputed as a matched run-level magnitude:
`0.5 * (abs(ON_run1 - OFF_run1) + abs(ON_run2 - OFF_run2))`.

These two subject-level quantities are compared with a paired one-sample test
on `run-level medication-effect magnitude - run-to-run baseline`, separately for
intra-ROI FC, between-ROI FC, and the primary intra-minus-between contrast.

Subjects included in the comparison: 17.
Connectivity metric: `pearson_fisher_z`.
Voxel selection: `weighted-vigour`.
