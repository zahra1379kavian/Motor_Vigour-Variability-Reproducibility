# Full vs Task-Only ROI Quantification

## Inputs
- Vigour network: `data/voxel_weights_task1_bold0.6_beta0.6_smooth1.25_gamma1.5_bold_thr90.html` selected overlay, with brainstem voxels suppressed.
- Task activation map: `data/z_valu_standard_glm.nii.gz` thresholded at z >= 3.5, with brainstem voxels suppressed.
- White matter excluded: Harvard-Oxford subcortical: Left/Right Cerebral White Matter.
- Atlas mode: AAL3v2 + Harvard-Oxford + residual nearest-label fill.

## Totals
- Vigour network: 5,957 non-white-matter voxels.
- Task activation map: 21,610 non-white-matter voxels.
- Same-voxel overlap: 503 voxels.
- Vigour-only: 5,454 voxels.
- Task-only: 21,107 voxels.
- Union: 27,064 voxels.

## Atlas Coverage
- AAL3v2 + Harvard-Oxford + residual nearest-label fill: union 27,064/27,064 (100.0%), vigour 100.0%, task 100.0%
- AAL3v2 coarse: union 25,810/27,064 (95.4%), vigour 98.1%, task 94.7%
- Harvard-Oxford cortical+subcortical: union 22,688/27,064 (83.8%), vigour 92.0%, task 80.9%
- Schaefer2018 100 parcels: union 20,670/27,064 (76.4%), vigour 83.9%, task 73.4%

## Region Sets
These sets are voxel-level summaries. A single ROI can contain both vigour-only and task-only voxels.
- regions_with_vigour_only_voxels: 23 regions - Parietal; Cerebellum; Temporal; Occipital; Frontal; Postcentral; Fusiform; Insula; Cingulate; Putamen; Precentral; Orbitofrontal; Supp_Motor_Area; ParaHippocampal; Paracentral_Lobule; Rolandic_Oper; Pallidum; Caudate; Hippocampus; Amygdala; Cerebral_Cortex; Thalamus; Olfactory
- regions_with_task_only_voxels: 20 regions - Parietal; Cerebellum; Temporal; Occipital; Frontal; Postcentral; Fusiform; Insula; Cingulate; Putamen; Precentral; Orbitofrontal; Supp_Motor_Area; ParaHippocampal; Paracentral_Lobule; Rolandic_Oper; Pallidum; Hippocampus; Cerebral_Cortex; Thalamus
- regions_with_same_voxel_overlap: 10 regions - Parietal; Cerebellum; Temporal; Occipital; Postcentral; Fusiform; Putamen; Precentral; Supp_Motor_Area; Paracentral_Lobule
- regions_with_both_maps_but_no_same_voxel_overlap: 10 regions - Frontal; Insula; Cingulate; Orbitofrontal; ParaHippocampal; Rolandic_Oper; Pallidum; Hippocampus; Cerebral_Cortex; Thalamus
- regions_present_only_in_vigour_map: 3 regions - Caudate; Amygdala; Olfactory
- regions_present_only_in_task_activation_map: 0 regions - None

## Highest Vigour-Only Percent of ROI
- Amygdala: 20.9% of ROI
- ParaHippocampal: 20.9% of ROI
- Caudate: 10.0% of ROI
- Orbitofrontal: 9.9% of ROI
- Supp_Motor_Area: 8.3% of ROI
- Paracentral_Lobule: 7.9% of ROI
- Olfactory: 7.8% of ROI
- Temporal: 7.4% of ROI

## Highest Task-Only Percent of ROI
- Putamen: 41.6% of ROI
- Pallidum: 32.7% of ROI
- Insula: 25.7% of ROI
- Postcentral: 24.2% of ROI
- Parietal: 23.9% of ROI
- Fusiform: 21.8% of ROI
- Occipital: 19.1% of ROI
- Cerebellum: 17.7% of ROI

## Highest Shared-Voxel Percent of ROI
- Putamen: 6.0% of ROI
- Paracentral_Lobule: 1.2% of ROI
- Cerebellum: 1.0% of ROI
- Occipital: 0.4% of ROI
- Postcentral: 0.3% of ROI
- Parietal: 0.3% of ROI
- Precentral: 0.3% of ROI
- Fusiform: 0.2% of ROI

## Outputs
- `figures/ablation/ablation_full_vs_task_only_roi_quantification_by_roi.csv`
- `figures/ablation/ablation_full_vs_task_only_roi_quantification_region_sets.csv`
- `figures/ablation/ablation_full_vs_task_only_roi_quantification_clean_table.csv`
- `figures/ablation/ablation_full_vs_task_only_roi_quantification_summary.png`
- `figures/ablation/ablation_full_vs_task_only_roi_quantification_clean_table.png`
- `figures/ablation/ablation_full_vs_task_only_roi_quantification_atlas_coverage.csv`
- `figures/ablation/ablation_full_vs_task_only_roi_quantification_metadata.json`