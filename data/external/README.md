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

