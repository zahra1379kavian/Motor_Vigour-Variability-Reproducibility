# Full vs Task-Only ROI Quantification

## Inputs
- Vigour network: `data/voxel_weights_task1_bold0.6_beta0.6_smooth1.25_gamma1.5_bold_thr90.html` selected overlay, with brainstem voxels suppressed.
- Task activation map: `data/z_valu_standard_glm.nii.gz` thresholded at z >= 3.5, with brainstem voxels suppressed.
- White matter excluded: Harvard-Oxford subcortical: Left/Right Cerebral White Matter.
- Atlas mode: AAL3v2 bilateral subregions + Harvard-Oxford exact-label fill + residual nearest-label fill.

## Totals
- Vigour network: 5,957 non-white-matter voxels.
- Task activation map: 21,610 non-white-matter voxels.
- Same-voxel overlap: 503 voxels.
- Vigour-only: 5,454 voxels.
- Task-only: 21,107 voxels.
- Union: 27,064 voxels.

## Atlas Coverage
- AAL3v2 bilateral subregions + Harvard-Oxford exact-label fill + residual nearest-label fill: union 27,064/27,064 (100.0%), vigour 100.0%, task 100.0%
- AAL3v2 coarse: union 25,810/27,064 (95.4%), vigour 98.1%, task 94.7%
- Harvard-Oxford cortical+subcortical: union 22,688/27,064 (83.8%), vigour 92.0%, task 80.9%
- Schaefer2018 100 parcels: union 20,670/27,064 (76.4%), vigour 83.9%, task 73.4%

## Region Sets
These sets are voxel-level summaries. A single ROI can contain both vigour-only and task-only voxels.
- regions_with_vigour_only_voxels: 59 regions - Postcentral; Occipital_Mid; Temporal_Inf; Cerebellum_6; SupraMarginal; Cerebellum_8; Temporal_Mid; Fusiform; Parietal_Sup; Precuneus; Parietal_Inf; Insula; Frontal_Mid_2; Putamen; Precentral; Cingulate_Mid; Cerebellum_4_5; Temporal_Sup; Frontal_Sup_2; Cerebellum_Crus1; Occipital_Sup; Supp_Motor_Area; Temporal_Pole_Mid; ParaHippocampal; Paracentral_Lobule; Frontal_Inf_Oper; Rolandic_Oper; Angular; Cuneus; Rectus; Pallidum; OFCpost; Cerebellum_9; OFCmed; Caudate; Hippocampus; Cerebellum_7b; Amygdala; Frontal_Sup_Medial; Cerebral_Cortex; Temporal_Pole_Sup; Temporal_Pole; Olfactory; Thal_VL; Postcentral_Gyrus; Temporal_Fusiform_Cortex_anterior_division; Precentral_Gyrus; Lateral_Occipital_Cortex_superior_division; Inferior_Temporal_Gyrus_posterior_division; ACC_sub; OFCant; Parahippocampal_Gyrus_anterior_division; Supramarginal_Gyrus_posterior_division; Thalamus; Frontal_Orbital_Cortex; Frontal_Med_Orb; Angular_Gyrus; Middle_Temporal_Gyrus_posterior_division; Subcallosal_Cortex
- regions_with_task_only_voxels: 87 regions - Postcentral; Occipital_Mid; Temporal_Inf; Cerebellum_6; SupraMarginal; Cerebellum_8; Temporal_Mid; Fusiform; Parietal_Sup; Precuneus; Parietal_Inf; Insula; Frontal_Mid_2; Putamen; Occipital_Inf; Precentral; Cingulate_Mid; Cerebellum_4_5; Lingual; Temporal_Sup; Frontal_Sup_2; Cerebellum_Crus1; Occipital_Sup; Supp_Motor_Area; ParaHippocampal; Paracentral_Lobule; Frontal_Inf_Oper; Rolandic_Oper; Angular; Calcarine; Frontal_Inf_Tri; Frontal_Inf_Orb_2; Vermis_4_5; Cuneus; Pallidum; OFCpost; Cerebellum_9; ACC_sup; Insular_Cortex; Vermis_8; Vermis_6; Vermis_9; Hippocampus; Cerebellum_7b; Frontal_Sup_Medial; Cerebral_Cortex; Temporal_Pole_Sup; Lateral_Occipital_Cortex_inferior_division; Heschl; Cerebellum_10; Temporal_Pole; Thal_VL; Postcentral_Gyrus; Cingulate_Gyrus_posterior_division; Superior_Parietal_Lobule; Vermis_7; Precentral_Gyrus; Supramarginal_Gyrus_anterior_division; Vermis_10; Lateral_Occipital_Cortex_superior_division; Cerebellum_Crus2; Thal_VPL; Thal_LGN; ACC_pre; OFCant; Supramarginal_Gyrus_posterior_division; Thalamus; Frontal_Orbital_Cortex; Lingual_Gyrus; Supracalcarine_Cortex; OFClat; Occipital_Fusiform_Gyrus; Precuneous_Cortex; Central_Opercular_Cortex; Cerebellum_3; Middle_Temporal_Gyrus_temporooccipital_part; Occipital_Pole; Thal_PuI; Thal_VA; Cingulate_Gyrus_anterior_division; Thal_AV; Frontal_Opercular_Cortex; Inferior_Temporal_Gyrus_temporooccipital_part; Parahippocampal_Gyrus_posterior_division; Parietal_Opercular_Cortex; Superior_Frontal_Gyrus; Thal_PuL
- regions_with_same_voxel_overlap: 19 regions - Postcentral; Occipital_Mid; Temporal_Inf; Cerebellum_6; Cerebellum_8; Fusiform; Parietal_Sup; Precuneus; Putamen; Precentral; Cerebellum_4_5; Temporal_Sup; Cerebellum_Crus1; Occipital_Sup; Supp_Motor_Area; Paracentral_Lobule; Cerebellum_9; Cerebellum_7b; Postcentral_Gyrus
- regions_with_both_maps_but_no_same_voxel_overlap: 26 regions - SupraMarginal; Temporal_Mid; Parietal_Inf; Insula; Frontal_Mid_2; Cingulate_Mid; Frontal_Sup_2; ParaHippocampal; Frontal_Inf_Oper; Rolandic_Oper; Angular; Cuneus; Pallidum; OFCpost; Hippocampus; Frontal_Sup_Medial; Cerebral_Cortex; Temporal_Pole_Sup; Temporal_Pole; Thal_VL; Precentral_Gyrus; Lateral_Occipital_Cortex_superior_division; OFCant; Supramarginal_Gyrus_posterior_division; Thalamus; Frontal_Orbital_Cortex
- regions_present_only_in_vigour_map: 14 regions - Temporal_Pole_Mid; Rectus; OFCmed; Caudate; Amygdala; Olfactory; Temporal_Fusiform_Cortex_anterior_division; Inferior_Temporal_Gyrus_posterior_division; ACC_sub; Parahippocampal_Gyrus_anterior_division; Frontal_Med_Orb; Angular_Gyrus; Middle_Temporal_Gyrus_posterior_division; Subcallosal_Cortex
- regions_present_only_in_task_activation_map: 42 regions - Occipital_Inf; Lingual; Calcarine; Frontal_Inf_Tri; Frontal_Inf_Orb_2; Vermis_4_5; ACC_sup; Insular_Cortex; Vermis_8; Vermis_6; Vermis_9; Lateral_Occipital_Cortex_inferior_division; Heschl; Cerebellum_10; Cingulate_Gyrus_posterior_division; Superior_Parietal_Lobule; Vermis_7; Supramarginal_Gyrus_anterior_division; Vermis_10; Cerebellum_Crus2; Thal_VPL; Thal_LGN; ACC_pre; Lingual_Gyrus; Supracalcarine_Cortex; OFClat; Occipital_Fusiform_Gyrus; Precuneous_Cortex; Central_Opercular_Cortex; Cerebellum_3; Middle_Temporal_Gyrus_temporooccipital_part; Occipital_Pole; Thal_PuI; Thal_VA; Cingulate_Gyrus_anterior_division; Thal_AV; Frontal_Opercular_Cortex; Inferior_Temporal_Gyrus_temporooccipital_part; Parahippocampal_Gyrus_posterior_division; Parietal_Opercular_Cortex; Superior_Frontal_Gyrus; Thal_PuL

## Highest Vigour-Only Percent of ROI
- ParaHippocampal: 23.5% of ROI
- Amygdala: 20.9% of ROI
- Temporal_Pole_Mid: 18.0% of ROI
- Temporal_Inf: 16.2% of ROI
- Rectus: 15.5% of ROI
- OFCmed: 14.5% of ROI
- OFCpost: 13.6% of ROI
- Temporal_Fusiform_Cortex_anterior_division: 11.8% of ROI

## Highest Task-Only Percent of ROI
- Cerebellum_3: 75.0% of ROI
- Thal_LGN: 58.8% of ROI
- Vermis_10: 53.6% of ROI
- Vermis_8: 50.8% of ROI
- Thal_PuI: 50.0% of ROI
- SupraMarginal: 49.7% of ROI
- Vermis_9: 48.0% of ROI
- Cerebellum_4_5: 45.9% of ROI

## Highest Shared-Voxel Percent of ROI
- Putamen: 6.0% of ROI
- Cerebellum_6: 2.4% of ROI
- Cerebellum_8: 2.3% of ROI
- Occipital_Mid: 1.5% of ROI
- Paracentral_Lobule: 1.2% of ROI
- Parietal_Sup: 1.0% of ROI
- Postcentral_Gyrus: 0.6% of ROI
- Cerebellum_4_5: 0.5% of ROI

## Outputs
- `figures/ablation/ablation_full_vs_task_only_roi_quantification_subregions_by_roi.csv`
- `figures/ablation/ablation_full_vs_task_only_roi_quantification_subregions_region_sets.csv`
- `figures/ablation/ablation_full_vs_task_only_roi_quantification_subregions_clean_table.csv`
- `figures/ablation/ablation_full_vs_task_only_roi_quantification_subregions_summary.png`
- `figures/ablation/ablation_full_vs_task_only_roi_quantification_subregions_clean_table.png`
- `figures/ablation/ablation_full_vs_task_only_roi_quantification_subregions_atlas_coverage.csv`
- `figures/ablation/ablation_full_vs_task_only_roi_quantification_subregions_metadata.json`