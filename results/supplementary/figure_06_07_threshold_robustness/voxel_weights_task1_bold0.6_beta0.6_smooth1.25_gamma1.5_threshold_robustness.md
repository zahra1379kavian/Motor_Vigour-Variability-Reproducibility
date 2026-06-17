# Threshold-Robustness Analysis

Input map: `data/voxel_weights_task1_bold0.6_beta0.6_smooth1.25_gamma1.5.nii.gz`

Reference visualization: p90-style thresholding of the final voxel-weight network.

## Method

- Thresholded nonzero voxel weights at percentiles: p80, p85, p90, p95, p97p5.
- Used p90 as the reference threshold and treated p85/p95 as the main relaxed/tightened sensitivity range.
- The displayed p90 montage from `data/voxel_weights_task1_bold0.6_beta0.6_smooth1.25_gamma1.5_bold_thr90.html` contains 6,797 voxels; this is shown separately from the raw percentile sweep.
- Assigned suprathreshold voxels to the current AAL3 atlas using AAL3v2 (AAL3v1.nii).
- Merged left/right AAL labels and repeated subparcels into coarser bilateral anatomical groups; the atlas map itself was not changed.
- Counted a group as reportable when it contained at least 25 suprathreshold voxels.

## Result

The main network is robust at the anatomical-group level across p85-p95: most p90 groups remain reportable after tightening to p95, while relaxing to p85 mainly expands already-present nodes.

At p90, the map contains 9,189 suprathreshold voxels, 625 connected components, and 19 reportable AAL groups. The atlas assigns 93.5% of p90 suprathreshold voxels; the remaining 595 voxels are kept in the CSV as unassigned but excluded from group-stability counts. Relaxing to p85 gives 20 reportable groups; tightening to p95 retains 84.2% of the p90 groups.

Stable p85-p95 groups: Amygdala, Caudate, Cerebellum, Frontal, Fusiform, Hippocampus, Occipital, Orbitofrontal, ParaHippocampal, Paracentral_Lobule, Parietal, Postcentral, Precentral, Putamen, Supp_Motor_Area, Temporal

Added when relaxed to p85: N_Acc

Dropped below report threshold when tightened to p95: Cingulate, Olfactory, Thalamus

Largest p90 groups by voxel count: Temporal (2200), Parietal (1142), Frontal (929), Cerebellum (659), Orbitofrontal (646), ParaHippocampal (448), Precentral (351), Fusiform (348), Supp_Motor_Area (318), Occipital (293), Postcentral (247), Caudate (217); plus 7 more

## Suggested Reporting

Report the p90 network as the primary visualization and include the robustness heatmap/table as a supplement. In text, emphasize anatomical-group stability rather than raw voxel overlap, because percentile thresholds are nested by construction. A clear phrasing is:

"Threshold sensitivity was assessed by repeating the network definition at p85, p90, and p95 of the nonzero voxel-weight distribution. The main p90 network was stable at the bilateral anatomical-group level: the same core AAL3-derived groups persisted across p85-p95, while threshold relaxation mainly expanded the network and threshold tightening removed smaller peripheral groups."

Then list the stable core and the threshold-sensitive nodes from the bullets above.

## Outputs

- `figures/voxel_weights_task1_bold0.6_beta0.6_smooth1.25_gamma1.5_threshold_robustness.png`
- `figures/voxel_weights_task1_bold0.6_beta0.6_smooth1.25_gamma1.5_threshold_robustness.pdf`
- `figures/voxel_weights_task1_bold0.6_beta0.6_smooth1.25_gamma1.5_threshold_robustness_atlas_regions.png`
- `figures/voxel_weights_task1_bold0.6_beta0.6_smooth1.25_gamma1.5_threshold_robustness_atlas_regions.pdf`
- `figures/voxel_weights_task1_bold0.6_beta0.6_smooth1.25_gamma1.5_threshold_robustness_summary.csv`
- `figures/voxel_weights_task1_bold0.6_beta0.6_smooth1.25_gamma1.5_threshold_robustness_regions.csv`
- `figures/voxel_weights_task1_bold0.6_beta0.6_smooth1.25_gamma1.5_threshold_robustness.json`
