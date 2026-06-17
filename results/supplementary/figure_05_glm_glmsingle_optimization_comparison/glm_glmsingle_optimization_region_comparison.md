# GLM, GLMsingle, and Optimization Region Comparison

## Inputs

- Standard GLM: `data/z_valu_standard_glm.nii.gz`
- GLMsingle Type A: `data/z_value_typaA.nii.gz`
- GLMsingle Type D: `data/z_value_typeD.nii.gz`
- Optimization weights: `data/voxel_weights_task1_bold0.6_beta0.6_smooth1.25_gamma1.5.nii.gz`

## Method

- All images were checked on the same grid and summarized in MNI 2 mm space.
- The analysis mask was defined as the union of nonzero GLMsingle Type A, GLMsingle Type D, and optimization-weight voxels (91,874 voxels). This avoids treating the standard-GLM full field of view as highlighted signal.
- Positive z-map highlights used Standard GLM: z >= 3.1; GLMsingle Type A: z >= 1.96; GLMsingle Type D: z >= 1.96. Negative and absolute-z summaries were also exported for sensitivity checks.
- The optimization map was summarized as a weight-importance map using p80, p90 nonzero-weight percentile thresholds, plus voxel-count-matched top-weight masks for each positive-z and absolute-z map.
- Atlas summaries used AAL3v2 (AAL3v1.nii), with the same bilateral coarse AAL3 grouping as the threshold-robustness figure.

## Primary Positive-Z vs Weight Summary

### Standard GLM positive z
- Highlighted voxels: 22,882
- Atlas-assigned fraction: 92.4%
- Reportable regions: 16
- Top regions: Cerebellum: 5,883 voxels (25.7% of highlighted map); Occipital: 3,372 voxels (14.7% of highlighted map); Parietal: 3,101 voxels (13.6% of highlighted map); Temporal: 2,022 voxels (8.8% of highlighted map); Frontal: 1,307 voxels (5.7% of highlighted map); Postcentral: 1,272 voxels (5.6% of highlighted map); Fusiform: 1,025 voxels (4.5% of highlighted map); Cingulate: 930 voxels (4.1% of highlighted map)
### GLMsingle Type A positive z
- Highlighted voxels: 26,123
- Atlas-assigned fraction: 92.4%
- Reportable regions: 16
- Top regions: Cerebellum: 5,399 voxels (20.7% of highlighted map); Occipital: 4,098 voxels (15.7% of highlighted map); Parietal: 3,147 voxels (12.0% of highlighted map); Precentral: 1,763 voxels (6.7% of highlighted map); Frontal: 1,763 voxels (6.7% of highlighted map); Temporal: 1,604 voxels (6.1% of highlighted map); Postcentral: 1,450 voxels (5.6% of highlighted map); Supp_Motor_Area: 1,104 voxels (4.2% of highlighted map)
### GLMsingle Type D positive z
- Highlighted voxels: 12,421
- Atlas-assigned fraction: 94.0%
- Reportable regions: 16
- Top regions: Occipital: 2,384 voxels (19.2% of highlighted map); Cerebellum: 2,048 voxels (16.5% of highlighted map); Parietal: 1,859 voxels (15.0% of highlighted map); Temporal: 964 voxels (7.8% of highlighted map); Frontal: 845 voxels (6.8% of highlighted map); Precentral: 834 voxels (6.7% of highlighted map); Postcentral: 748 voxels (6.0% of highlighted map); Supp_Motor_Area: 468 voxels (3.8% of highlighted map)
### Optimization p80
- Highlighted voxels: 18,375
- Atlas-assigned fraction: 93.2%
- Reportable regions: 22
- Top regions: Temporal: 3,662 voxels (19.9% of highlighted map); Parietal: 2,449 voxels (13.3% of highlighted map); Cerebellum: 2,008 voxels (10.9% of highlighted map); Frontal: 1,711 voxels (9.3% of highlighted map); Occipital: 1,063 voxels (5.8% of highlighted map); Orbitofrontal: 1,063 voxels (5.8% of highlighted map); Fusiform: 739 voxels (4.0% of highlighted map); ParaHippocampal: 587 voxels (3.2% of highlighted map)
### Optimization p90
- Highlighted voxels: 9,189
- Atlas-assigned fraction: 93.5%
- Reportable regions: 19
- Top regions: Temporal: 2,200 voxels (23.9% of highlighted map); Parietal: 1,142 voxels (12.4% of highlighted map); Frontal: 929 voxels (10.1% of highlighted map); Cerebellum: 659 voxels (7.2% of highlighted map); Orbitofrontal: 646 voxels (7.0% of highlighted map); ParaHippocampal: 448 voxels (4.9% of highlighted map); Precentral: 351 voxels (3.8% of highlighted map); Fusiform: 348 voxels (3.8% of highlighted map)

## Primary Pairwise Spatial Overlap

- Optimization p80 vs Optimization p90: Dice=0.667, Jaccard=0.500, shared=9,189; Temporal (2200); Parietal (1142); Frontal (929); Cerebellum (659); Orbitofrontal (646); ParaHippocampal (448)
- Standard GLM positive z vs GLMsingle Type A positive z: Dice=0.639, Jaccard=0.470, shared=15,660; Cerebellum (4382); Occipital (2453); Parietal (2260); Temporal (1365); Postcentral (834); Insula (642)
- GLMsingle Type A positive z vs GLMsingle Type D positive z: Dice=0.556, Jaccard=0.385, shared=10,715; Occipital (2224); Cerebellum (1702); Parietal (1648); Temporal (813); Precentral (766); Frontal (659)
- Standard GLM positive z vs GLMsingle Type D positive z: Dice=0.433, Jaccard=0.277, shared=7,651; Cerebellum (1658); Occipital (1540); Parietal (1368); Temporal (796); Postcentral (410); Frontal (345)
- GLMsingle Type A positive z vs Optimization p80: Dice=0.207, Jaccard=0.115, shared=4,595; Cerebellum (915); Parietal (838); Occipital (530); Precentral (441); Postcentral (361); Temporal (282)
- Standard GLM positive z vs Optimization p80: Dice=0.182, Jaccard=0.100, shared=3,754; Cerebellum (1069); Parietal (724); Occipital (471); Temporal (374); Postcentral (244); Fusiform (225)
- GLMsingle Type D positive z vs Optimization p80: Dice=0.165, Jaccard=0.090, shared=2,540; Parietal (555); Cerebellum (417); Occipital (341); Precentral (260); Postcentral (190); Temporal (188)
- GLMsingle Type A positive z vs Optimization p90: Dice=0.106, Jaccard=0.056, shared=1,870; Parietal (371); Cerebellum (349); Precentral (261); Postcentral (184); Occipital (139); Supp_Motor_Area (138)
- GLMsingle Type D positive z vs Optimization p90: Dice=0.100, Jaccard=0.052, shared=1,077; Parietal (255); Cerebellum (163); Precentral (157); Occipital (100); Postcentral (99); Supp_Motor_Area (86)
- Standard GLM positive z vs Optimization p90: Dice=0.086, Jaccard=0.045, shared=1,383; Cerebellum (416); Parietal (250); Occipital (127); Putamen (119); Postcentral (116); Temporal (105)

## Outputs

- `figures/glm_glmsingle_optimization_region_comparison_summary.csv`
- `figures/glm_glmsingle_optimization_region_comparison_regions.csv`
- `figures/glm_glmsingle_optimization_region_comparison_overlaps.csv`
- `figures/glm_glmsingle_optimization_region_comparison_region_by_method.csv`
- `figures/glm_glmsingle_optimization_region_comparison_region_heatmap.png`
- `figures/glm_glmsingle_optimization_region_comparison_region_heatmap.pdf`
- `figures/glm_glmsingle_optimization_region_comparison_overlap_heatmap.png`
- `figures/glm_glmsingle_optimization_region_comparison_overlap_heatmap.pdf`
- `figures/glm_glmsingle_optimization_region_comparison.json`
