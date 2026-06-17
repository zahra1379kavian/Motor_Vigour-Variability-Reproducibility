# Task-Activation Intra-ROI vs Between-ROI FC Run-Baseline Comparison

This companion analysis uses the task-activation ROI masks from
`med_effects_task_activation.py`: standard-GLM z-map voxels thresholded at
z >= 3.1, unit voxel weights, and the same
lateralized AAL grouping used for the task-activation medication-change figure.

For each subject/session/run, intra-ROI voxel-pair FC and between-ROI FC were
computed separately.

The plotted run-to-run baseline is a within-session variability magnitude:
`0.5 * (abs(OFF_run2 - OFF_run1) + abs(ON_run2 - ON_run1))`.

The plotted medication effect is a matched run-level magnitude:
`0.5 * (abs(ON_run1 - OFF_run1) + abs(ON_run2 - OFF_run2))`.

The paired contrast tests `run-level medication-effect magnitude -
run-to-run baseline` across subjects, separately for intra-ROI FC,
between-ROI FC, and the intra-minus-between contrast.

Subjects included in the comparison: 17.
Connectivity metric: `pearson_fisher_z`.
