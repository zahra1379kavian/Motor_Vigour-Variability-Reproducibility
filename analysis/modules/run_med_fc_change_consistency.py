#!/usr/bin/env python3
"""Does medication change a subject's FC network, and is the change consistent?

Two questions, two tests, on the per-subject OFF (ses-1) / ON (ses-2) networks:

Q1 -- DOES MEDICATION CHANGE THE NETWORK (beyond run-to-run noise)?
  We need a within-state noise floor. By default this script reads a
  run-labelled beta-series table and estimates FC separately for run 1 and run
  2 within each session:
      within-state distance  = d(OFF_run1, OFF_run2) and d(ON_run1, ON_run2)
      between-state distance = d(OFF_run*, ON_run*)
  If medication changed nothing, ON and OFF differ only by the same noise that
  separates two runs of one session, so between == within.  A subject whose
  between-state distance exceeds their within-state distance has an FC network
  that changed with medication beyond run-to-run variability.  We test
  (between - within) across subjects (paired t + sign-flip perm).
  Distance = correlation distance (1 - Pearson r of the edge vectors).

  A legacy --split-half mode is retained for the previous random half-split
  analysis of already-concatenated per-session ROI time series.

Q2 -- IS THE CHANGE CONSISTENT ACROSS SUBJECTS?
  Each subject's full-session change vector is delta_i = FC_ON_i - FC_OFF_i.
  If subjects reorganise their network the same way, their delta vectors point
  in a common direction.  We quantify this as the mean pairwise Pearson
  correlation between subjects' delta vectors, and as a leave-one-out
  consistency (corr of each subject's delta with the mean of the others).
  Null: randomly sign-flip each subject's delta vector (preserves each
  subject's change magnitude, randomises direction) -> distribution of the
  consistency statistic under "no shared direction".

Outputs (under med_fc_change_consistency/):
  * change_vs_noise_subject.csv     - per subject within/between distance
  * change_consistency_subject.csv  - per subject leave-one-out consistency
  * change_consistency_summary.csv  - group stats per FC metric
  * med_fc_change_consistency.png/.pdf
  * med_fc_change_consistency_method.md
"""


import argparse
import importlib.util
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats as sstats
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
# Default = vigour-weighted network. Overridable via CLI to e.g. the
# task-activation ROI definition in figures/med_effects_task_activation.
TS_DIR = ROOT / "results" / "main" / "figure_06a_medication_vigour_network" / "roi_timeseries"
OUT = ROOT / "results" / "supplementary" / "figure_09_medication_fc_consistency"
DEFAULT_TRIAL_TABLE = (
    ROOT / "data" / "processed" / "gvs_connectivity" / "vigour_network" / "per_trial_roi_betas.csv"
)
GVS_METRIC_SCRIPT = (
    ROOT / "analysis" / "gvs_connectivity_coactivation" / "run_edge_connectivity_metric_sensitivity.py"
)

OFF_SESSION = "ses-1"
ON_SESSION = "ses-2"
RUN_1 = "1"
RUN_2 = "2"
N_SPLITS = 50          # random half-splits for the noise floor
N_PERM = 10000
TRIAL_META_COLUMNS = {
    "subject", "session", "medication", "run", "condition_code",
    "condition_label", "trial_in_condition",
}


def _load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Could not import {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


gvs = _load(GVS_METRIC_SCRIPT, "gvs_edge_metric_sensitivity")

FC_METRICS = [
    ("pearson_r", gvs.pearson_edges),
    ("spearman_rho", gvs.spearman_edges),
    ("partial_corr_ledoitwolf", gvs.partial_corr_edges),
    ("mutual_info_quantile", gvs.mutual_information_edges),
]


def load_timeseries(subject, session):
    path = TS_DIR / f"{subject}_{session}.csv"
    if not path.exists():
        return None
    df = pd.read_csv(path)
    df = df.loc[~df.isna().all(axis=1)]
    return df.to_numpy(dtype=float)


def _session_label(value):
    text = str(value).strip()
    if text.startswith("ses-"):
        return text
    if text.endswith(".0"):
        text = text[:-2]
    return f"ses-{text}"


def _run_label(value):
    text = str(value).strip()
    if text.endswith(".0"):
        text = text[:-2]
    return text


def load_run_timeseries(trial_table):
    """Load run-level ROI beta series from a long table with subject/session/run columns."""
    if not trial_table.exists():
        raise FileNotFoundError(f"Missing run-labelled beta-series table: {trial_table}")
    df = pd.read_csv(trial_table)
    required = {"subject", "session", "run"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"{trial_table} is missing columns: {', '.join(sorted(missing))}")
    roi_cols = [c for c in df.columns if c not in TRIAL_META_COLUMNS]
    if len(roi_cols) < 2:
        raise ValueError(f"{trial_table} has fewer than two ROI beta-series columns")

    work = df.copy()
    work["_session_label"] = work["session"].map(_session_label)
    work["_run_label"] = work["run"].map(_run_label)
    run_ts = {}
    for (subject, session, run), group in work.groupby(
        ["subject", "_session_label", "_run_label"], sort=True
    ):
        roi_df = group.loc[:, roi_cols].apply(pd.to_numeric, errors="coerce")
        roi_df = roi_df.loc[~roi_df.isna().all(axis=1)]
        if roi_df.empty:
            continue
        run_ts[(str(subject), str(session), str(run))] = roi_df.to_numpy(dtype=float)
    return roi_cols, run_ts


def complete_subjects():
    subs = sorted({p.name.split("_ses")[0] for p in TS_DIR.glob("*.csv")})
    return [
        s for s in subs
        if (TS_DIR / f"{s}_{OFF_SESSION}.csv").exists()
        and (TS_DIR / f"{s}_{ON_SESSION}.csv").exists()
    ]


def complete_run_subjects(run_ts):
    subjects = sorted({subject for subject, _, _ in run_ts})
    required = {
        (OFF_SESSION, RUN_1), (OFF_SESSION, RUN_2),
        (ON_SESSION, RUN_1), (ON_SESSION, RUN_2),
    }
    keep = []
    for subject in subjects:
        available = {(session, run) for s, session, run in run_ts if s == subject}
        if required <= available:
            keep.append(subject)
    return keep


def corr_distance(a, b):
    mask = np.isfinite(a) & np.isfinite(b)
    x, y = a[mask], b[mask]
    if x.size < 3 or np.std(x) == 0 or np.std(y) == 0:
        return np.nan
    return 1.0 - float(np.corrcoef(x, y)[0, 1])


def run_distances(run_ts, subject, fc_fn):
    """Return (within_state_dist, between_state_dist) from run 1/run 2 networks."""
    off_edges = {
        run: fc_fn(run_ts[(subject, OFF_SESSION, run)])
        for run in (RUN_1, RUN_2)
    }
    on_edges = {
        run: fc_fn(run_ts[(subject, ON_SESSION, run)])
        for run in (RUN_1, RUN_2)
    }
    within = [
        corr_distance(off_edges[RUN_1], off_edges[RUN_2]),
        corr_distance(on_edges[RUN_1], on_edges[RUN_2]),
    ]
    between = [
        corr_distance(off_edges[off_run], on_edges[on_run])
        for off_run in (RUN_1, RUN_2)
        for on_run in (RUN_1, RUN_2)
    ]
    return float(np.nanmean(within)), float(np.nanmean(between))


def split_distances(off_ts, on_ts, fc_fn, rng):
    """Return (within_state_dist, between_state_dist) averaged over splits."""
    within, between = [], []
    for _ in range(N_SPLITS):
        # independent random half-split of each session's time points
        def halves(ts):
            n = ts.shape[0]
            idx = rng.permutation(n)
            h = n // 2
            return ts[idx[:h]], ts[idx[h:2 * h]]

        o1, o2 = halves(off_ts)
        n1, n2 = halves(on_ts)
        fo1, fo2 = fc_fn(o1), fc_fn(o2)
        fn1, fn2 = fc_fn(n1), fc_fn(n2)
        within.append(corr_distance(fo1, fo2))
        within.append(corr_distance(fn1, fn2))
        between.extend([
            corr_distance(fo1, fn1), corr_distance(fo1, fn2),
            corr_distance(fo2, fn1), corr_distance(fo2, fn2),
        ])
    return float(np.nanmean(within)), float(np.nanmean(between))


def concatenate_session_runs(run_ts, subject, session):
    return np.vstack([
        run_ts[(subject, session, RUN_1)],
        run_ts[(subject, session, RUN_2)],
    ])


def paired_signflip_p(delta, rng):
    delta = delta[np.isfinite(delta)]
    n = delta.size
    if n < 3 or np.std(delta) == 0:
        return np.nan
    obs = abs(np.mean(delta))
    signs = rng.choice([-1.0, 1.0], size=(N_PERM, n))
    null = np.abs((signs * delta[None, :]).mean(axis=1))
    return float((np.sum(null >= obs - 1e-12) + 1) / (N_PERM + 1))


def consistency_stats(deltas, rng):
    """deltas: subjects x edges (ON - OFF). Mean pairwise corr + LOO + perm."""
    # drop edges that are NaN for any subject
    good = np.isfinite(deltas).all(axis=0)
    D = deltas[:, good]
    n = D.shape[0]
    # column-center? No -- we want directional alignment, use raw Pearson per pair.
    R = np.corrcoef(D)
    iu = np.triu_indices(n, k=1)
    mean_pairwise = float(np.nanmean(R[iu]))
    # leave-one-out: corr(subject_i delta, mean of others)
    loo = []
    for i in range(n):
        others = np.delete(np.arange(n), i)
        m = D[others].mean(axis=0)
        x = D[i]
        if np.std(x) > 0 and np.std(m) > 0:
            loo.append(float(np.corrcoef(x, m)[0, 1]))
    loo = np.array(loo)
    mean_loo = float(np.mean(loo))
    # sign-flip null on mean pairwise corr
    null = np.empty(N_PERM)
    for k in range(N_PERM):
        s = rng.choice([-1.0, 1.0], size=(n, 1))
        Dp = D * s
        Rp = np.corrcoef(Dp)
        null[k] = np.nanmean(Rp[iu])
    p_perm = float((np.sum(null >= mean_pairwise - 1e-12) + 1) / (N_PERM + 1))
    # parametric test that LOO consistency > 0
    if loo.size > 2 and np.std(loo) > 0:
        t_loo, p_loo = sstats.ttest_1samp(loo, 0.0)
        p_loo_greater = p_loo / 2 if t_loo > 0 else 1 - p_loo / 2
    else:
        t_loo, p_loo_greater = np.nan, np.nan
    return dict(mean_pairwise_corr=mean_pairwise, p_perm_consistency=p_perm,
                mean_loo_corr=mean_loo, t_loo=float(t_loo),
                p_loo_greater=float(p_loo_greater), loo=loo)


def main():
    global TS_DIR, OUT
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--trial-table", type=Path, default=DEFAULT_TRIAL_TABLE,
                        help="run-labelled ROI beta-series table; default uses the vigour-weighted table")
    parser.add_argument("--ts-dir", type=Path, default=TS_DIR,
                        help="directory of per-subject/session ROI timeseries CSVs")
    parser.add_argument("--out", type=Path, default=OUT,
                        help="output directory")
    parser.add_argument("--split-half", action="store_true",
                        help="legacy mode: use random split halves of concatenated session time series")
    parser.add_argument("--merge-panel", nargs=2, action="append", default=[],
                        metavar=("LABEL", "CONSISTENCY_DIR"),
                        help="make a merged figure from existing output dirs; repeat for each row")
    parser.add_argument("--merged-out", type=Path, default=OUT,
                        help="directory for merged-network figure when --merge-panel is used")
    args = parser.parse_args()
    TS_DIR = args.ts_dir
    OUT = args.out

    if args.merge_panel:
        _make_merged_figure([(label, Path(path)) for label, path in args.merge_panel],
                            args.merged_out)
        return

    source = "split_half" if args.split_half else "run1_vs_run2"
    print(f"OUT={OUT}")
    if args.split_half:
        print(f"TS_DIR={TS_DIR}\nNoise floor: random split halves")
    else:
        print(f"TRIAL_TABLE={args.trial_table}\nNoise floor: run 1 vs run 2")

    OUT.mkdir(parents=True, exist_ok=True)
    if args.split_half:
        subjects = complete_subjects()
        run_ts = None
    else:
        _, run_ts = load_run_timeseries(args.trial_table)
        subjects = complete_run_subjects(run_ts)
    n_sub = len(subjects)
    print(f"Subjects (OFF+ON): {n_sub}")

    noise_rows, loo_rows, summary_rows = [], [], []

    for fc_name, fc_fn in FC_METRICS:
        if args.split_half:
            off = {s: load_timeseries(s, OFF_SESSION) for s in subjects}
            on = {s: load_timeseries(s, ON_SESSION) for s in subjects}
        else:
            off = {s: concatenate_session_runs(run_ts, s, OFF_SESSION) for s in subjects}
            on = {s: concatenate_session_runs(run_ts, s, ON_SESSION) for s in subjects}

        # ---- Q1: change beyond noise ----
        rng = np.random.default_rng(0)
        within_d, between_d = [], []
        for s in subjects:
            if args.split_half:
                w, b = split_distances(off[s], on[s], fc_fn, rng)
            else:
                w, b = run_distances(run_ts, s, fc_fn)
            within_d.append(w)
            between_d.append(b)
            noise_rows.append(dict(fc_metric=fc_name, subject=s,
                                   noise_floor=source,
                                   within_state_dist=w, between_state_dist=b,
                                   excess=b - w))
        within_d = np.array(within_d)
        between_d = np.array(between_d)
        excess = between_d - within_d
        t_e, p_e = sstats.ttest_1samp(excess[np.isfinite(excess)], 0.0)
        p_e_greater = p_e / 2 if t_e > 0 else 1 - p_e / 2
        d_e = float(np.nanmean(excess) / np.nanstd(excess, ddof=1))
        p_e_perm = paired_signflip_p(excess, np.random.default_rng(7))
        n_pos = int(np.sum(excess > 0))

        # ---- Q2: consistency of the change direction ----
        deltas = np.vstack([fc_fn(on[s]) - fc_fn(off[s]) for s in subjects])
        cs = consistency_stats(deltas, np.random.default_rng(11))
        for s, lo in zip(subjects, cs["loo"]):
            loo_rows.append(dict(fc_metric=fc_name, subject=s, loo_corr=lo))

        summary_rows.append(dict(
            fc_metric=fc_name, n_subjects=n_sub, noise_floor=source,
            mean_within_dist=float(np.nanmean(within_d)),
            mean_between_dist=float(np.nanmean(between_d)),
            mean_excess=float(np.nanmean(excess)), excess_cohens_d=d_e,
            excess_t=float(t_e), excess_p_greater=float(p_e_greater),
            excess_p_signflip=p_e_perm, n_subjects_changed=n_pos,
            mean_pairwise_delta_corr=cs["mean_pairwise_corr"],
            consistency_p_perm=cs["p_perm_consistency"],
            mean_loo_corr=cs["mean_loo_corr"],
            loo_t=cs["t_loo"], loo_p_greater=cs["p_loo_greater"]))

        print(f"  {fc_name:24s} | Q1 within={np.nanmean(within_d):.3f} "
              f"between={np.nanmean(between_d):.3f} excess d={d_e:.2f} "
              f"p={p_e_greater:.3g} perm={p_e_perm:.3g} changed={n_pos}/{n_sub} "
              f"|| Q2 pairwise r={cs['mean_pairwise_corr']:.3f} "
              f"p={cs['p_perm_consistency']:.3g} LOO r={cs['mean_loo_corr']:.3f}")

    noise = pd.DataFrame(noise_rows)
    loo = pd.DataFrame(loo_rows)
    summary = pd.DataFrame(summary_rows)
    noise.to_csv(OUT / "change_vs_noise_subject.csv", index=False)
    loo.to_csv(OUT / "change_consistency_subject.csv", index=False)
    summary.to_csv(OUT / "change_consistency_summary.csv", index=False)

    _make_figure(noise, summary)
    (OUT / "med_fc_change_consistency_method.md").write_text(
        _method_note(n_sub, source=source, trial_table=args.trial_table if not args.split_half else None)
    )
    print("\nSummary:")
    cols = ["fc_metric", "mean_within_dist", "mean_between_dist", "mean_excess",
            "excess_cohens_d", "excess_p_greater", "excess_p_signflip",
            "n_subjects_changed", "mean_pairwise_delta_corr",
            "consistency_p_perm", "mean_loo_corr", "loo_p_greater"]
    print(summary[cols].to_string(index=False))


def _noise_style(noise_floor):
    if noise_floor == "run1_vs_run2":
        return "within-state run 1 vs 2", ["run 1\nvs run 2", "ON vs\nOFF"], "run-to-run noise"
    return "within-state split halves", ["within\nhalves", "ON vs\nOFF"], "split-half noise"


def _summary_noise_floor(summary):
    if "noise_floor" in summary.columns and not summary["noise_floor"].dropna().empty:
        return str(summary["noise_floor"].dropna().iloc[0])
    return "split_half"


def _load_edge_summary(panel_dir):
    path = panel_dir.parent / "within_subject_med_fc" / "within_subject_metric_summary.csv"
    if not path.exists():
        return None
    return pd.read_csv(path)


def _edge_title_line(edge_summary, fc, g):
    if edge_summary is None or "metric" not in edge_summary.columns:
        return (f"consistency r={g['mean_pairwise_delta_corr']:.3f} "
                f"p={g['consistency_p_perm']:.3g}")
    row = edge_summary.loc[edge_summary["metric"] == fc]
    if row.empty:
        return (f"consistency r={g['mean_pairwise_delta_corr']:.3f} "
                f"p={g['consistency_p_perm']:.3g}")
    row = row.iloc[0]
    return (f"edge-FDR q<0.05: t={int(row['n_sig_edges_fdr_ttest'])}, "
            f"perm={int(row['n_sig_edges_fdr_perm'])}")


def _draw_metric_axis(ax, noise, summary, fc, *, show_ylabel, show_legend, edge_summary=None):
    d = noise[noise.fc_metric == fc]
    g = summary[summary.fc_metric == fc].iloc[0]
    noise_floor = str(g["noise_floor"]) if "noise_floor" in g.index else _summary_noise_floor(summary)
    noise_label, ticklabels, _ = _noise_style(noise_floor)
    w, b = d["within_state_dist"].to_numpy(), d["between_state_dist"].to_numpy()
    for wi, bi in zip(w, b):
        ax.plot([0, 1], [wi, bi], color="0.7", lw=0.7, zorder=1)
    ax.scatter(np.zeros_like(w), w, color="#3182bd", s=16, zorder=2,
               label=noise_label)
    ax.scatter(np.ones_like(b), b, color="#e6550d", s=16, zorder=2,
               label="ON vs OFF")
    ax.set_xticks([0, 1])
    ax.set_xticklabels(ticklabels, fontsize=8)
    star = "*" if (np.isfinite(g["excess_p_signflip"]) and g["excess_p_signflip"] < 0.05) else ""
    ax.set_title(f"{fc}\nchange>noise: d={g['excess_cohens_d']:.2f} "
                 f"p={g['excess_p_signflip']:.3g}{star}\n"
                 f"{_edge_title_line(edge_summary, fc, g)}", fontsize=7.5)
    if show_ylabel:
        ax.set_ylabel("correlation distance", fontsize=8)
    ax.tick_params(labelsize=7)
    if show_legend:
        ax.legend(fontsize=6, loc="upper left")


def _make_figure(noise, summary):
    fcs = list(dict.fromkeys(summary["fc_metric"]))
    fig, axes = plt.subplots(1, len(fcs), figsize=(3.2 * len(fcs), 4.0), squeeze=False)
    _, _, noise_title = _noise_style(_summary_noise_floor(summary))
    edge_summary = _load_edge_summary(OUT)
    for c, fc in enumerate(fcs):
        ax = axes[0][c]
        _draw_metric_axis(ax, noise, summary, fc, show_ylabel=(c == 0),
                          show_legend=(c == 0), edge_summary=edge_summary)
    fig.suptitle(f"Does medication change the FC network (vs {noise_title}), and which "
                 f"edges change consistently? (n={int(summary['n_subjects'].max())})",
                 fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    fig.savefig(OUT / "med_fc_change_consistency.png", dpi=200)
    fig.savefig(OUT / "med_fc_change_consistency.pdf")
    plt.close(fig)


def _make_merged_figure(panel_specs, out_dir):
    if not panel_specs:
        raise ValueError("No panels were provided for merged figure")
    panels = []
    for label, panel_dir in panel_specs:
        noise_path = panel_dir / "change_vs_noise_subject.csv"
        summary_path = panel_dir / "change_consistency_summary.csv"
        if not noise_path.exists() or not summary_path.exists():
            raise FileNotFoundError(f"Missing consistency CSVs under {panel_dir}")
        panels.append((label, pd.read_csv(noise_path), pd.read_csv(summary_path),
                       _load_edge_summary(panel_dir)))

    fcs = list(dict.fromkeys(panels[0][2]["fc_metric"]))
    fig, axes = plt.subplots(len(panels), len(fcs), figsize=(3.2 * len(fcs), 3.25 * len(panels)),
                             squeeze=False)
    for r, (label, noise, summary, edge_summary) in enumerate(panels):
        for c, fc in enumerate(fcs):
            _draw_metric_axis(axes[r][c], noise, summary, fc,
                              show_ylabel=(c == 0), show_legend=(c == 0),
                              edge_summary=edge_summary)
        axes[r][0].text(-0.24, 1.28, label, transform=axes[r][0].transAxes,
                        ha="left", va="bottom", fontsize=13, fontweight="bold")
    fig.tight_layout(h_pad=3.0)
    out_dir.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_dir / "med_fc_change_consistency_merged_networks.png", dpi=200)
    fig.savefig(out_dir / "med_fc_change_consistency_merged_networks.pdf")
    plt.close(fig)


def _method_note(n_subjects, *, source, trial_table):
    fc = "\n".join(f"  - `{n}`" for n, _ in FC_METRICS)
    if source == "run1_vs_run2":
        q1 = f"""**Q1 -- change beyond run-to-run noise.** Run-labelled ROI beta series were
read from `{trial_table}`. FC was estimated separately for run 1 and run 2
within each OFF and ON session. The within-state run distance is the mean of
d(OFF_run1,OFF_run2) and d(ON_run1,ON_run2). The between-state distance is the
mean of the four OFF-run versus ON-run distances. Distance = correlation
distance (1 - Pearson r of edge vectors). Per subject we test
(between - within) > 0 across subjects with a paired t-test and a sign-flip
permutation test."""
    else:
        q1 = f"""**Q1 -- change beyond split-half noise.** Each session's time series is split
into two random halves ({N_SPLITS} splits) and FC is estimated from each half,
so every network estimate rests on the same number of time points. The
within-state distance d(OFF_h1,OFF_h2)/d(ON_h1,ON_h2) is the split-half noise
floor; the between-state distance d(OFF_h*,ON_h*) adds any true medication
effect. Distance = correlation distance (1 - Pearson r of edge vectors). Per
subject we average across splits, then test (between - within) > 0 across
subjects with a paired t-test and a sign-flip permutation test."""
    return f"""# Does medication change the FC network, and is it consistent?

Subjects with OFF (ses-1) and ON (ses-2) sessions: {n_subjects}.

{q1}

**Q2 -- consistency across subjects.** Each subject's change vector is
delta_i = FC_ON_i - FC_OFF_i. We report the mean pairwise Pearson correlation
between subjects' delta vectors and a leave-one-out consistency (corr of each
subject's delta with the average delta of the others). Null: randomly sign-flip
each subject's delta vector ({N_PERM} permutations), which keeps each subject's
change magnitude but destroys any shared direction.

FC definitions:
{fc}
"""


if __name__ == "__main__":
    main()
