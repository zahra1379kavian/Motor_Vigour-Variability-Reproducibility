# Does medication change the FC network, and is it consistent?

Subjects with OFF (ses-1) and ON (ses-2) sessions: 17.

**Q1 -- change beyond run-to-run noise.** Run-labelled ROI beta series were
read from `/home/zkavian/Master_thesis_final/Claude_results/group_analyses/analysis3_lme/per_trial_roi_betas.csv`. FC was estimated separately for run 1 and run 2
within each OFF and ON session. The within-state run distance is the mean of
d(OFF_run1,OFF_run2) and d(ON_run1,ON_run2). The between-state distance is the
mean of the four OFF-run versus ON-run distances. Distance = correlation
distance (1 - Pearson r of edge vectors). Per subject we test
(between - within) > 0 across subjects with a paired t-test and a sign-flip
permutation test.

**Q2 -- consistency across subjects.** Each subject's change vector is
delta_i = FC_ON_i - FC_OFF_i. We report the mean pairwise Pearson correlation
between subjects' delta vectors and a leave-one-out consistency (corr of each
subject's delta with the average delta of the others). Null: randomly sign-flip
each subject's delta vector (10000 permutations), which keeps each subject's
change magnitude but destroys any shared direction.

FC definitions:
  - `pearson_r`
  - `spearman_rho`
  - `partial_corr_ledoitwolf`
  - `mutual_info_quantile`
