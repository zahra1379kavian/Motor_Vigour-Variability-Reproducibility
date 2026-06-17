# Ablation Study Analysis

## Inputs

- Main map: `data/voxel_weights_task1_bold0.6_beta0.6_smooth1.25_gamma1.5.nii.gz`; summarized as voxels above the 90th percentile of nonzero weights (top 10%).
- Ablation maps: `data/ablation`; existing `*_bold_thr90.nii.gz` maps were treated as already thresholded.
- Correlation/generalization metrics: `data/ablation/slurm-11460024.out` and `data/ablation/slurm-11445550.out`.

## Outputs

- `ablation_fold_metrics.csv`: fold-level correlations, losses, term contributions, and balanced scores.
- `ablation_metric_summary.csv`: candidate-level means and SEMs.
- `ablation_spatial_overlap.csv`: Dice/Jaccard/center-of-mass distances relative to the main full p90 map.
- `ablation_map_summary.csv`: selected-voxel counts and cluster summaries.
- `ablation_roi_regions.csv`: AAL coarse ROI composition by map.
- `ablation_publication_summary.{png,pdf}`: compact multi-panel publication figure.
- `ablation_map_montage.{png,pdf}`: axial slice overview of available maps.
- `ablation_full_vs_task_only_anatomy.{png,pdf}`: Figure-3-style sagittal/coronal/axial contour montage for the supplied full-model and task-only thresholded HTML maps.
- `ablation_full_vs_task_only_anatomy_summary.csv`: voxel counts and overlap for the focused full-vs-task-only anatomy figure.

## Notes

- The final-weight SLURM log contains the `task=1, bold=0.6, beta=0.6, smooth=1.25, gamma=1.5` full model plus no-task/no-BOLD/no-beta and single-term baselines, but no final-weight no-smooth metrics.
- `slurm-11445550.out` is mixed; only the parameter-independent task-only and no-objective baselines were retained. The unrelated `task=1, bold=1, beta=0.75, smooth=1.8` sweep was excluded.
- Balanced scores were computed within the final-weight analysis group using inverse ranges of candidate-mean score components.
- The focused full-vs-task-only anatomy figure uses the selected-voxel overlay embedded in the full-model thresholded HTML map and `data/z_valu_standard_glm.nii.gz` thresholded at z >= 3.5; brainstem contours are suppressed, filled blue overlays show full-model-only voxels, filled vermillion overlays show standard-GLM-only voxels, and filled green overlays with a white halo show overlap with stronger line weight over motor ROIs.

## Balanced-Score Weights

| analysis_group   |   corr_weight |   corr_gap_weight |   loss_gap_weight | definition                                                             |
|:-----------------|--------------:|------------------:|------------------:|:-----------------------------------------------------------------------|
| final_0p6        |       24.5821 |           33.9328 |       3.31479e-06 | inverse range of candidate-mean score components within analysis_group |

## Missing Expected Maps

| analysis_group   | candidate_label   | candidate_id                           | expected_path                                                                                               | exists   | note                                                                                        |
|:-----------------|:------------------|:---------------------------------------|:------------------------------------------------------------------------------------------------------------|:---------|:--------------------------------------------------------------------------------------------|
| final_0p6        | No smoothness     | task1_bold0.6_beta0.6_smooth0_gamma1.5 | data/ablation/voxel_weights_mean_foldavg_sub9_ses1_task1_bold0.6_beta0.6_smooth0_gamma1.5_bold_thr90.nii.gz | False    | Final no-smooth map and metrics were not available in the retained final-parameter results. |

## Metric Summary

| analysis_group   | candidate_label        |   abs_test_corr_mean |   corr_generalization_gap_mean |   relative_loss_gap_mean |   balanced_score_mean |
|:-----------------|:-----------------------|---------------------:|-------------------------------:|-------------------------:|----------------------:|
| final_0p6        | No task                |              0.03771 |                        0.03367 |                0.103749  |              0.215527 |
| final_0p6        | BOLD-only              |              0.06026 |                        0.05393 |                0.205501  |              0.34868  |
| final_0p6        | Beta-only              |              0.07839 |                        0.06314 |                0.0245529 |              0.215527 |
| final_0p6        | No objective penalties |              0.07499 |                        0.04782 |                0.475942  |             -0.220743 |
| final_0p6        | Smooth-only            |              0.07636 |                        0.04722 |                0         |             -0.274782 |
| final_0p6        | Full model             |              0.0664  |                        0.05934 |                0.271468  |              0.381322 |
| final_0p6        | No beta stability      |              0.06836 |                        0.0601  |                0.412915  |              0.358931 |
| final_0p6        | No BOLD stability      |              0.05808 |                        0.04841 |                0.229821  |              0.21496  |
| final_0p6        | Task-only reference    |              0.05919 |                        0.05063 |           301678         |              1.263    |

## Spatial Summary

| label                  |   map_voxels |   shared_voxels |        dice |   center_of_mass_distance_mm |
|:-----------------------|-------------:|----------------:|------------:|-----------------------------:|
| No task                |         7200 |            2675 | 0.326438    |                     25.3969  |
| BOLD-only              |         5115 |            2000 | 0.279642    |                     32.0592  |
| Beta-only              |         6695 |              14 | 0.00176278  |                     33.5661  |
| No objective penalties |         5160 |             655 | 0.0912956   |                     10.0916  |
| Smooth-only            |         8517 |               2 | 0.000225912 |                     13.325   |
| Full ablation map      |         6252 |            6103 | 0.790493    |                      6.05978 |
| No beta stability      |         6154 |            5928 | 0.77273     |                      6.19358 |
| No BOLD stability      |         6104 |            4291 | 0.561172    |                     16.0547  |
| Task-only reference    |         4640 |            3473 | 0.502278    |                      9.71204 |
