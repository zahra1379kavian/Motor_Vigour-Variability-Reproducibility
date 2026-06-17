# Behavioural Supplementary Methods and Results

## Methods

Behavioural analyses used the first six columns of each revised metrics file: 1/PT, 1/RT, 1/MT, 1/(RT+MT), Vmax, and Pmax. The first four columns were reciprocal time measures and were converted to milliseconds for reporting. Vmax and Pmax were kept on their raw feature scales.

Each GVS condition contributed 20 trials per session. Trials 1-10 were assigned to run 1 and trials 11-20 to run 2. Actual block position within each run was read from `data/gvs_order_by_subject_session_run.tsv`.

Catch trials were identified from `res.catchtrials` and verified against reward code -5. Go trials were non-catch trials. Missed trials were go trials marked as no-squeeze in `res.sqrwd_nosq`. Valid RT trials were go trials with finite positive source 1/RT; invalid RT trials were non-missed go trials without finite positive 1/RT.

RT variability was defined as RMSSD in milliseconds, computed from adjacent valid RT pairs within subject-session-run-GVS blocks. This excludes transitions across GVS blocks, rest gaps, and run boundaries.

## Results

After excluding PSPD017 session 1, the analysis retained 5368 valid RT trials from 35 sessions and 18 subjects. Mean RT was 380.5 ms (SD across trials 101.1 ms).

RT was retained as the primary behavioural vigour measure because it is a direct, interpretable response-initiation latency with low missingness and a stable mapping to the trial-wise neural analyses. In the feature comparison table, RT missing rate was 0.053; run1-run2 r=0.81 (n=312). Force features were retained as manipulation checks rather than the primary vigour outcome.

Medication changed mean RT by ON-OFF=-13.6 ms (95% CI -27.5, 0.3; p = 0.054; n=17). Medication changed RT RMSSD by ON-OFF=-2.1 ms (95% CI -18.2, 14.0; p = 0.787; n=17).

GVS effects on RT were small. The smallest FDR-adjusted condition-level RT effect was Beta with q=0.624 and p = 0.078.

Across blocks, mean RT RMSSD was 83.8 ms, mean lag-1 autocorrelation was -0.147, and mean linear drift was -0.23 ms/trial.

Reward-residualized neural sensitivity (mean_projection_delta_vs_reward_residualized_mean_rt_delta) gave r=0.071, p = 0.238, n=280.
Reward-residualized neural sensitivity (projection_variability_delta_vs_reward_residualized_rt_rmssd_delta) gave r=0.106, p = 0.076, n=280.

The generated CSV tables contain the full trial flow, feature comparison, model coefficients, GVS condition tests, temporal structure, and subject-level heterogeneity summaries.
