# External Data Requirements

Some analyses require subject-level inputs that are not packaged in this
repository.

Expected locations used by the copied scripts:

```text
/mnt/TeamShare/Data_Masterfile/H20-00572_All-Dressed/Zahra-Thesis-Data/fmri_opt_group/results_beta_preprocessed
/mnt/TeamShare/Data_Masterfile/H20-00572_All-Dressed/Zahra-Thesis-Data/fmri_opt_group/behaviour
/mnt/TeamShare/Data_Masterfile/H20-00572_All-Dressed/Zahra-Thesis-Data/fmri_opt_group/data/group_concat
/usr/local/fsl/data/atlases/Cerebellum/Cerebellum-MNIfnirt-maxprob-thr25-2mm.nii.gz
```

Nilearn atlases may also be cached under:

```text
/home/zkavian/nilearn_data
```

If these paths differ on another machine, pass the corresponding command-line
arguments shown by `python scripts/<script_name>.py --help`.

Supplementary Figures 13 and 17 also use large working-analysis arrays that are
not packaged here. To rerun the copied scripts with default arguments, place the
following files under `data/external/` or pass replacement paths:

```text
GVS_projection_BOLD/
Task_map_BOLD_trials.npy
active_bold_group.npy
active_flat_indices__group.npy
concat_manifest_group.tsv
z_valu_standard_glm_zgte3.1_active_bold.npy
ablation/voxel_weights_mean_foldavg_sub9_ses1_task0_bold0_beta0_smooth0_gamma1.5_bold_thr90_active_bold.npy
ablation/voxel_weights_mean_foldavg_sub9_ses1_task1_bold0_beta0.6_smooth1.25_gamma1.5_bold_thr90_active_bold.npy
ablation/voxel_weights_mean_foldavg_sub9_ses1_task1_bold0.6_beta0_smooth1.25_gamma1.5_bold_thr90_active_bold.npy
ablation/voxel_weights_mean_foldavg_sub9_ses1_task1_bold0.6_beta0.6_smooth1.25_gamma1.5_bold_thr90_active_bold.npy
```
