# Behaviour Supplement Assumptions and Definitions

## Source Files

- Behaviour metrics: `/mnt/TeamShare/Data_Masterfile/H20-00572_All-Dressed/PRECISIONSTIM_PD_Data_Results/Behaviour/Behaviour_metrics_revised`
- Consolidated metadata: `/mnt/TeamShare/Data_Masterfile/H20-00572_All-Dressed/PRECISIONSTIM_PD_Data_Results/Behaviour/Consolidated_behav_data`
- GVS order inventory: `/home/zkavian/Master_thesis_final/data/gvs_order_by_subject_session_run.tsv`
- GVS parameter metadata: `/mnt/TeamShare/Data_Masterfile/H20-00572_All-Dressed/AllDressed_WorkOnData/Sepideh/PredictiveModel/gvs_params.csv`

The Windows source path supplied by the user maps to `/mnt/TeamShare/Data_Masterfile/.../PRECISIONSTIM_PD_Data_Results/Behaviour/Behaviour_metrics_revised` in this environment.

## Trial Definitions

- Total trials: all 20 entries per GVS condition split into run 1 trials 1-10 and run 2 trials 11-20.
- Catch trials: trials with `res.catchtrials == 1`; this was cross-checked against reward code -5.
- Go trials: all non-catch trials.
- Missed trials: go trials marked as no squeeze in `res.sqrwd_nosq == 1`.
- Invalid RT trials: go trials that were not missed but did not have finite positive source `1/RT`.
- Valid RT trials: go trials with finite positive source `1/RT`; RT is reported as `1000 / (1/RT)` in ms.

Very-late responses with finite positive RT were retained as valid RT trials to match the existing behavioural scripts, which use finite positive RT rather than reward-success filtering.

## Exclusions

- PSPD017 session 1 is marked and excluded from group analyses because one recording run was excessively noisy.
- The primary analysis uses the session numbering in the GVS/fMRI inventory: OFF=session 1 and ON=session 2.
- The older `figures/reward_effects` output used a subject-specific PSPD017 medication/session override, so reward effects are recomputed here from raw files for consistency.

## RT Variability

RT variability is root mean squared successive difference (RMSSD) in milliseconds. It is computed only for adjacent valid RT pairs within the same subject, session, run, and GVS block. Transitions across GVS blocks, rest periods, and run boundaries are excluded.

## Non-Dominant Task Hand

The prompt notes that one right-handed participant used the left hand because dominant-hand tremor was severe. The local behavioural metrics, consolidated metadata, GVS order files, and available demographics tables did not identify which subject this was. The script therefore exposes `--non-dominant-task-hand-subject PSPD###` so the exact participant can be marked when the ID is known.
