# Intra-ROI vs Between-ROI FC Method

For each subject/session, beta-series were extracted from task-activation voxels defined by the standard-GLM z map thresholded at z >= 3.1. Because the task map is binary after thresholding and has no optimization weights, ROI beta-series are unweighted means of selected voxel beta values.

Intra-ROI FC was defined as the mean Pearson correlation between voxel beta-series within the same ROI. Voxel-pair correlations were Fisher z transformed before averaging, with equal weight for each voxel pair. Between-ROI FC was computed from the unweighted mean ROI beta-series in the same Fisher-z Pearson-correlation scale. Medication effects were evaluated within complete subjects as ON minus OFF, and the primary comparison was (ON - OFF intra-ROI FC) - (ON - OFF between-ROI FC).
