# Supplementary Figure 4: Standard GLM FEAT Result

This folder contains the compact provenance package for Supplementary Figure 4.

Source figure:

```text
/home/zkavian/fsl_glm/outputs/feat/mixed_model.gfeat/cope1.feat/rendered_thresh_zstat1_z3p1_brain_slices_whitebg_top_RL.png
```

Curated figure export:

```text
figures/supplementary/supp_figure_04_standard_glm_feat_brain_slices.png
figures/supplementary/supp_figure_04_standard_glm_feat_brain_slices.pdf
```

Included provenance:

- `feat_design/`: mixed-model and cope-level FEAT design files.
- `feat_stats/`: compact thresholded/statistical NIfTI outputs and cluster/local-maxima tables.
- `feat_reports/`: FEAT HTML reports.
- `provenance/`: mixed-model FSF, run table, trial-selection summary, EV files, and generated first-level/session fixed FSFs.
- `analysis/fsl_glm/`: Python code copied from the source FSL GLM workflow.

Large FEAT volumes such as `filtered_func_data.nii.gz`, `res4d.nii.gz`, and
`weights1.nii.gz` were intentionally not copied.
