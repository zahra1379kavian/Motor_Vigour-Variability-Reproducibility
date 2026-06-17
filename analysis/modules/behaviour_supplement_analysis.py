#!/usr/bin/env python3
"""Generate behavioural supplementary tables, figures, and text.

This script is intentionally self-contained so the Supplementary Material
analysis can be rerun from one command. It reads the revised behavioural
metrics and consolidated metadata, applies the documented subject/session
exclusion, and writes all derived outputs into figures/Behav_supp.
"""


import argparse
import json
import math
import os
import re
import shutil
import warnings
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", "/tmp/matplotlib-behav-supp")

import matplotlib

warnings.filterwarnings("ignore", message="Unable to import Axes3D.*")
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.text import Text
import numpy as np
import pandas as pd
from scipy import stats
from scipy.io import loadmat

try:
    import statsmodels.formula.api as smf
except ImportError:  # pragma: no cover
    smf = None


ROOT = Path(__file__).resolve().parents[2]
OUT_DIR = ROOT / "results" / "supplementary" / "figure_01_02_behavior"
DEFAULT_SOURCE_ROOT = Path("/mnt/TeamShare/Data_Masterfile/H20-00572_All-Dressed/" "AllDressed_WorkOnData/Sepideh")
DEFAULT_BEHAVIOUR_ROOT = Path("/mnt/TeamShare/Data_Masterfile/H20-00572_All-Dressed/" "PRECISIONSTIM_PD_Data_Results/Behaviour")
DEFAULT_METRICS_DIR = DEFAULT_BEHAVIOUR_ROOT / "Behaviour_metrics_revised"
DEFAULT_CONSOLIDATED_DIR = DEFAULT_BEHAVIOUR_ROOT / "Consolidated_behav_data"
DEFAULT_GVS_ORDER = ROOT / "data" / "metadata" / "gvs_order_by_subject_session_run.tsv"
DEFAULT_GVS_PARAMS = DEFAULT_SOURCE_ROOT / "PredictiveModel" / "gvs_params.csv"
DEFAULT_DOMINANT_HAND = Path("/mnt/TeamShare/Data_Masterfile/H20-00572_All-Dressed/" "Zahra-Thesis-Data/fsl_gl-temporary/analysis_qc/" "mixed_model_responder_summary_with_dominant_hand.tsv")
DEFAULT_PRIOR_GVS_DIR = ROOT / "results" / "supplementary" / "figure_11_gvs_projection_rt"
DEFAULT_PRIOR_REWARD_DIR = ROOT / "results" / "supplementary" / "figure_03_reward_effects"
DEFAULT_NEURAL_MEAN = ROOT / "data" / "external" / "B_projected_signal_raw.csv"
DEFAULT_NEURAL_VARIABILITY = (ROOT / "data" / "external" / "B_projected_signal_consecutive_diff_subject_session.csv")

SESSION_BY_MEDICATION = {"OFF": 1, "ON": 2}
EXCLUDED_SUBJECT_SESSIONS = {("PSPD017", 1): "excluded: session 1 had one excessively noisy recording run"}
LOW_REWARD_CODES = {0, 1}
HIGH_REWARD_CODES = {5}
CATCH_REWARD_CODE = -5
REST_AFTER_BLOCKS = {3, 6}

PAPER_STYLE = {"font.family": "sans-serif", "font.sans-serif": ["Liberation Sans", "Arial", "DejaVu Sans"], "pdf.fonttype": 42, "ps.fonttype": 42, "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "grid.color": "#E5E7EB", "grid.linewidth": 0.7, "grid.alpha": 0.8}

BOLD_PAPER_STYLE = {**PAPER_STYLE, "font.size": 16, "font.weight": "bold", "axes.labelsize": 18, "axes.labelweight": "bold", "axes.titlesize": 18, "axes.titleweight": "bold", "xtick.labelsize": 15, "ytick.labelsize": 15, "legend.fontsize": 16}


def _bold_figure_text(fig):
    for text in fig.findobj(match=Text):
        text.set_fontweight("bold")


class Measure:
    def __init__(self, name, source_column, source_label, output_column, unit, kind, axis_label):
        self.name = name
        self.source_column = source_column
        self.source_label = source_label
        self.output_column = output_column
        self.unit = unit
        self.kind = kind
        self.axis_label = axis_label


MEASURES = [
    Measure("PT", 0, "1/PT", "PT_ms", "ms", "inverse_time", "PT (ms)"),
    Measure("RT", 1, "1/RT", "RT_ms", "ms", "inverse_time", "RT (ms)"),
    Measure("MT", 2, "1/MT", "MT_ms", "ms", "inverse_time", "MT (ms)"),
    Measure("RT+MT", 3, "1/(RT+MT)", "RT_plus_MT_ms", "ms", "inverse_time", "RT + MT (ms)"),
    Measure("Vmax", 4, "Vmax", "Vmax", "raw feature units", "raw", "Vmax (raw units)"),
    Measure("Pmax", 5, "Pmax", "Pmax", "raw feature units", "raw", "Pmax (raw units)"),
]

STIM_SHORT = {1: "Sham", 2: "Pink noise", 3: "DC+1", 4: "DC-1", 5: "Delta", 6: "Theta", 7: "Alpha", 8: "Beta", 9: "Gamma"}


def _parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--metrics-dir", type=Path, default=DEFAULT_METRICS_DIR)
    parser.add_argument("--consolidated-dir", type=Path, default=DEFAULT_CONSOLIDATED_DIR)
    parser.add_argument("--gvs-order", type=Path, default=DEFAULT_GVS_ORDER)
    parser.add_argument("--gvs-params", type=Path, default=DEFAULT_GVS_PARAMS)
    parser.add_argument("--dominant-hand-table", type=Path, default=DEFAULT_DOMINANT_HAND)
    parser.add_argument("--prior-gvs-dir", type=Path, default=DEFAULT_PRIOR_GVS_DIR)
    parser.add_argument("--prior-reward-dir", type=Path, default=DEFAULT_PRIOR_REWARD_DIR)
    parser.add_argument("--neural-mean-table", type=Path, default=DEFAULT_NEURAL_MEAN)
    parser.add_argument("--neural-variability-table", type=Path, default=DEFAULT_NEURAL_VARIABILITY)
    parser.add_argument("--out-dir", type=Path, default=OUT_DIR)
    parser.add_argument("--non-dominant-task-hand-subject", default=None, help=("Optional subject ID, for example PSPD001, if the right-handed " "participant who used the left hand is known outside the local files."))
    parser.add_argument("--skip-models", action="store_true", help="Build data/figures but skip model fitting.")
    return parser.parse_args()


def _format_p(value):
    if not np.isfinite(value):
        return "p = NA"
    if value < 0.001:
        return "p < 0.001"
    return f"p = {value:.3f}"


def _format_num(value, digits=2):
    if not np.isfinite(value):
        return "NA"
    return f"{value:.{digits}f}"


def _subject_number(subject):
    match = re.search(r"(\d+)", str(subject))
    if match is None:
        raise ValueError(f"Could not parse subject number from {subject!r}")
    return int(match.group(1))


def _pspd_id(subject):
    return f"PSPD{_subject_number(str(subject)):03d}"


def _subpd_id(subject):
    return f"sub-pd{_subject_number(str(subject)):03d}"


def _opposite_hand(hand):
    hand = "" if hand is None or (isinstance(hand, float) and np.isnan(hand)) else str(hand).upper()
    if hand == "R":
        return "L"
    if hand == "L":
        return "R"
    return ""


def _bh_fdr(p_values):
    p_values = np.asarray(p_values, dtype=float)
    out = np.full_like(p_values, np.nan, dtype=float)
    finite = np.isfinite(p_values)
    if not np.any(finite):
        return out
    p = p_values[finite]
    order = np.argsort(p)
    ranked = p[order]
    n = ranked.size
    adjusted = ranked * n / np.arange(1, n + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    restored = np.empty_like(adjusted)
    restored[order] = np.clip(adjusted, 0.0, 1.0)
    out[finite] = restored
    return out


def _safe_ttest_1samp(values):
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    row = {"n": int(values.size), "mean": float(np.mean(values)) if values.size else np.nan, "median": float(np.median(values)) if values.size else np.nan, "sd": float(np.std(values, ddof=1)) if values.size > 1 else np.nan, "sem": float(stats.sem(values)) if values.size > 1 else np.nan, "ci95_low": np.nan, "ci95_high": np.nan, "t": np.nan, "p": np.nan, "cohen_dz": np.nan}
    if values.size > 1:
        t_res = stats.ttest_1samp(values, 0.0)
        ci_low, ci_high = stats.t.interval(0.95, values.size - 1, loc=row["mean"], scale=row["sem"])
        row.update({"ci95_low": float(ci_low), "ci95_high": float(ci_high), "t": float(t_res.statistic), "p": float(t_res.pvalue), "cohen_dz": float(row["mean"] / row["sd"]) if row["sd"] and row["sd"] > 0 else np.nan})
    return row


def _safe_pearson(x, y):
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    keep = np.isfinite(x) & np.isfinite(y)
    if np.count_nonzero(keep) < 3:
        return {"n": int(np.count_nonzero(keep)), "r": np.nan, "p": np.nan}
    r, p = stats.pearsonr(x[keep], y[keep])
    return {"n": int(np.count_nonzero(keep)), "r": float(r), "p": float(p)}


def _mean_sem(values):
    arr = values.to_numpy(dtype=float)
    arr = arr[np.isfinite(arr)]
    if arr.size == 0:
        return np.nan, np.nan
    return float(np.mean(arr)), float(stats.sem(arr)) if arr.size > 1 else np.nan


def _mean_ci(values, confidence=0.95):
    arr = np.asarray(values, dtype=float)
    arr = arr[np.isfinite(arr)]
    if arr.size == 0:
        return np.nan, np.nan, np.nan
    mean = float(np.mean(arr))
    if arr.size < 2:
        return mean, np.nan, np.nan
    sem = float(stats.sem(arr))
    half_width = float(stats.t.ppf(0.5 + confidence / 2.0, arr.size - 1) * sem)
    return mean, mean - half_width, mean + half_width


def _save_figure(fig, stem, dpi=300):
    stem.parent.mkdir(parents=True, exist_ok=True)
    png = stem.with_suffix(".png")
    pdf = stem.with_suffix(".pdf")
    fig.savefig(png, dpi=dpi, bbox_inches="tight", pad_inches=0.04)
    fig.savefig(pdf, bbox_inches="tight", pad_inches=0.04)
    plt.close(fig)
    return [png, pdf]


def _write_csv(df, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)
    return path


def _write_markdown_table(df, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    try:
        text = df.to_markdown(index=False)
    except Exception:
        text = df.to_csv(index=False)
    path.write_text(text + "\n", encoding="utf-8")
    return path


def _load_gvs_params(path):
    if path.exists():
        params = pd.read_csv(path)
        params = params.loc[params["stim_id"].between(1, 9)].copy()
    else:
        params = pd.DataFrame({"stim_id": list(range(1, 10))})
    params["stim_id"] = params["stim_id"].astype(int)
    if "stim_label" not in params:
        params["stim_label"] = params["stim_id"].map(STIM_SHORT)
    for col in ["family", "main_freq_hz", "carrier_hz", "env_hz", "band", "phase_band", "amp_band"]:
        if col not in params:
            params[col] = np.nan
    params["gvs_condition"] = params["stim_id"].map(lambda x: f"gvs-{int(x):02d}")
    params["gvs_label"] = params["stim_id"].map(STIM_SHORT).fillna(params["stim_label"].astype(str))
    return params[["stim_id", "gvs_condition", "gvs_label", "stim_label", "family", "main_freq_hz", "carrier_hz", "env_hz", "band", "phase_band", "amp_band"]].copy()


def _load_dominant_hand(path, non_dominant_subject):
    if path.exists():
        hand = pd.read_csv(path, sep="\t")
        hand["subject"] = hand["subject"].map(_pspd_id)
        hand = hand.loc[:, ["subject", "dominant_hand"]].copy()
    else:
        hand = pd.DataFrame(columns=["subject", "dominant_hand"])
    if non_dominant_subject:
        non_dominant_subject = _pspd_id(non_dominant_subject)
    hand["task_hand"] = hand["dominant_hand"]
    hand["non_dominant_task_hand"] = False
    if non_dominant_subject:
        if non_dominant_subject not in set(hand["subject"]):
            hand = pd.concat([hand, pd.DataFrame([{"subject": non_dominant_subject, "dominant_hand": "", "task_hand": "L"}])], ignore_index=True)
        mask = hand["subject"].eq(non_dominant_subject)
        hand.loc[mask, "task_hand"] = hand.loc[mask, "dominant_hand"].map(_opposite_hand)
        hand.loc[mask & hand["task_hand"].eq(""), "task_hand"] = "L"
        hand.loc[mask, "non_dominant_task_hand"] = True
    hand["hand_note"] = np.where(hand["non_dominant_task_hand"], "right-handed participant used the left hand because dominant-hand tremor was severe", "")
    return hand


def _load_order_map(path):
    order_map = {}
    if not path.exists():
        return order_map
    order = pd.read_csv(path, sep="\t")
    order = order.loc[order.get("status", "ok").astype(str).eq("ok")].copy()
    for row in order.itertuples(index=False):
        subject = _pspd_id(row.subject_id)
        medication = str(row.session).upper()
        run = int(row.run)
        for position in range(1, 10):
            value = getattr(row, f"stim_order_{position}", np.nan)
            if np.isfinite(value):
                order_map[(subject, medication, run, int(value))] = position
    return order_map


def _load_metric_cells(path):
    mat = loadmat(path, squeeze_me=False, struct_as_record=False)
    if "behav_metrics" not in mat:
        raise ValueError(f"{path} does not contain behav_metrics")
    cells = np.asarray(mat["behav_metrics"], dtype=object).ravel()
    metrics = [np.asarray(cell, dtype=float) for cell in cells]
    if len(metrics) < 9:
        raise ValueError(f"{path} has {len(metrics)} metric cells; expected at least 9")
    return metrics[:9]


def _load_res_fields(path):
    mat = loadmat(path, squeeze_me=False, struct_as_record=False)
    if "res" not in mat:
        raise ValueError(f"{path} does not contain res")
    res = mat["res"][0, 0]
    fields = {}
    for name in ["catchtrials", "goodtrials", "man_rej", "others", "reward", "sqrwd_nosq", "sqrwd_prem", "sqrwd_success", "sqrwd_vlate", "sqtype"]:
        fields[name] = np.asarray(getattr(res, name), dtype=float) if hasattr(res, name) else np.full((20, 9), np.nan)
    return fields


def _source_to_value(value, measure):
    if not np.isfinite(value):
        return np.nan
    if measure.kind == "inverse_time":
        if value <= 0:
            return np.nan
        return 1000.0 / value
    return float(value)


def _valid_source_value(value, measure):
    if not np.isfinite(value):
        return False
    if measure.kind == "inverse_time":
        return value > 0
    return True


def _parse_metrics_name(path):
    match = re.match(r"^(PSPD\d+)_(OFF|ON)_behav_metrics\.mat$", path.name)
    if match is None:
        raise ValueError(f"Unrecognized metrics filename: {path.name}")
    return _pspd_id(match.group(1)), match.group(2)


def _combine_notes(values):
    cleaned = []
    for value in values:
        if isinstance(value, str) and value and value not in cleaned:
            cleaned.append(value)
    return "; ".join(cleaned)


def build_trial_table(metrics_dir, consolidated_dir, gvs_params, order_map, hand_table):
    params = gvs_params.set_index("stim_id").to_dict("index")
    hand = hand_table.set_index("subject").to_dict("index")
    rows = []
    metric_paths = sorted(metrics_dir.glob("PSPD*_behav_metrics.mat"))
    if not metric_paths:
        raise FileNotFoundError(f"No behavioural metrics files found under {metrics_dir}")

    for metric_path in metric_paths:
        subject, medication = _parse_metrics_name(metric_path)
        consolidated_path = consolidated_dir / metric_path.name.replace("_behav_metrics.mat", "_consolidated_behavdata.mat")
        if not consolidated_path.exists():
            raise FileNotFoundError(f"Missing consolidated metadata for {metric_path.name}: {consolidated_path}")

        session = SESSION_BY_MEDICATION[medication]
        excluded = (subject, session) in EXCLUDED_SUBJECT_SESSIONS
        exclusion_reason = EXCLUDED_SUBJECT_SESSIONS.get((subject, session), "")
        metric_cells = _load_metric_cells(metric_path)
        res = _load_res_fields(consolidated_path)
        hand_row = hand.get(subject, {})

        for stim_id, metric in enumerate(metric_cells, start=1):
            if metric.ndim != 2 or metric.shape[1] < 6:
                raise ValueError(f"{metric_path} stim {stim_id} has shape {metric.shape}; expected >=6 columns")
            stim_meta = params.get(stim_id, {})
            for trial_in_condition in range(1, min(metric.shape[0], 20) + 1):
                idx = trial_in_condition - 1
                run = 1 if trial_in_condition <= 10 else 2
                trial_in_run = ((trial_in_condition - 1) % 10) + 1
                block_position = order_map.get((subject, medication, run, stim_id), stim_id)
                run_trial_index = (block_position - 1) * 10 + trial_in_run
                reward_code = int(res["reward"][idx, stim_id - 1]) if np.isfinite(res["reward"][idx, stim_id - 1]) else np.nan
                catch_field = res["catchtrials"][idx, stim_id - 1]
                catch_trial = bool(catch_field == 1) if np.isfinite(catch_field) else reward_code == CATCH_REWARD_CODE
                go_trial = not catch_trial
                missed_trial = bool(go_trial and res["sqrwd_nosq"][idx, stim_id - 1] == 1)
                rt_source = float(metric[idx, 1])
                rt_source_valid = np.isfinite(rt_source) and rt_source > 0
                valid_rt = bool(go_trial and rt_source_valid)
                invalid_rt = bool(go_trial and (not missed_trial) and (not rt_source_valid))
                reward_level = "high" if reward_code in HIGH_REWARD_CODES else "low" if reward_code in LOW_REWARD_CODES else "catch"

                record = {
                    "subject": subject,
                    "subject_label": _subpd_id(subject),
                    "subject_number": _subject_number(subject),
                    "session": session,
                    "medication": medication,
                    "run": run,
                    "stim_id": stim_id,
                    "gvs_condition": stim_meta.get("gvs_condition", f"gvs-{stim_id:02d}"),
                    "gvs_label": stim_meta.get("gvs_label", STIM_SHORT.get(stim_id, f"GVS{stim_id}")),
                    "stim_label": stim_meta.get("stim_label", STIM_SHORT.get(stim_id, f"GVS{stim_id}")),
                    "family": stim_meta.get("family", np.nan),
                    "main_freq_hz": stim_meta.get("main_freq_hz", np.nan),
                    "carrier_hz": stim_meta.get("carrier_hz", np.nan),
                    "env_hz": stim_meta.get("env_hz", np.nan),
                    "band": stim_meta.get("band", np.nan),
                    "block_position": int(block_position),
                    "precedes_rest": bool(block_position in REST_AFTER_BLOCKS),
                    "follows_rest": bool(block_position in {block + 1 for block in REST_AFTER_BLOCKS}),
                    "trial_in_condition": int(trial_in_condition),
                    "trial_in_run": int(trial_in_run),
                    "run_trial_index": int(run_trial_index),
                    "reward_code": reward_code,
                    "reward_level": reward_level,
                    "reward_high": 1 if reward_level == "high" else 0 if reward_level == "low" else np.nan,
                    "catch_trial": catch_trial,
                    "go_trial": go_trial,
                    "missed_trial": missed_trial,
                    "invalid_rt_trial": invalid_rt,
                    "valid_rt_trial": valid_rt,
                    "goodtrial_flag": int(res["goodtrials"][idx, stim_id - 1])
                    if np.isfinite(res["goodtrials"][idx, stim_id - 1])
                    else np.nan,
                    "manual_reject_flag": int(res["man_rej"][idx, stim_id - 1])
                    if np.isfinite(res["man_rej"][idx, stim_id - 1])
                    else np.nan,
                    "premature_flag": int(res["sqrwd_prem"][idx, stim_id - 1])
                    if np.isfinite(res["sqrwd_prem"][idx, stim_id - 1])
                    else np.nan,
                    "very_late_flag": int(res["sqrwd_vlate"][idx, stim_id - 1])
                    if np.isfinite(res["sqrwd_vlate"][idx, stim_id - 1])
                    else np.nan,
                    "success_flag": int(res["sqrwd_success"][idx, stim_id - 1])
                    if np.isfinite(res["sqrwd_success"][idx, stim_id - 1])
                    else np.nan,
                    "sqtype": int(res["sqtype"][idx, stim_id - 1]) if np.isfinite(res["sqtype"][idx, stim_id - 1]) else np.nan,
                    "excluded_session": bool(excluded),
                    "exclusion_reason": exclusion_reason,
                    "dominant_hand": hand_row.get("dominant_hand", ""),
                    "task_hand": hand_row.get("task_hand", hand_row.get("dominant_hand", "")),
                    "non_dominant_task_hand": bool(hand_row.get("non_dominant_task_hand", False)),
                    "hand_note": hand_row.get("hand_note", ""),
                    "metrics_file": str(metric_path),
                    "consolidated_file": str(consolidated_path),
                }
                for measure in MEASURES:
                    source_value = float(metric[idx, measure.source_column])
                    record[f"source_{measure.source_label}"] = source_value
                    record[measure.output_column] = _source_to_value(source_value, measure)
                    record[f"valid_{measure.name.replace('+', '_plus_')}"] = bool(go_trial and _valid_source_value(source_value, measure))
                rows.append(record)

    trials = pd.DataFrame(rows)
    trials = trials.sort_values(["subject_number", "session", "run", "block_position", "trial_in_run", "stim_id"]).reset_index(drop=True)
    return trials


def build_trial_flow(trials):
    rows = []
    for key, group in trials.groupby(["subject", "session", "medication", "run"], sort=True):
        subject, session, medication, run = key
        notes = _combine_notes([str(group["exclusion_reason"].dropna().iloc[0]) if group["excluded_session"].any() else "", str(group["hand_note"].dropna().iloc[0]) if group["non_dominant_task_hand"].any() else ""])
        rows.append(
            {
                "Subject": subject,
                "Session": int(session),
                "Medication": medication,
                "Run": int(run),
                "Total trials": int(len(group)),
                "Go trials": int(group["go_trial"].sum()),
                "Catch trials": int(group["catch_trial"].sum()),
                "Missed trials": int(group["missed_trial"].sum()),
                "Invalid RT": int(group["invalid_rt_trial"].sum()),
                "Valid RT trials": int(group["valid_rt_trial"].sum()),
                "Session status": "excluded" if group["excluded_session"].any() else "included",
                "Notes": notes,
            }
        )
    return pd.DataFrame(rows)


def build_measure_long(trials):
    pieces = []
    base_cols = ["subject", "subject_label", "subject_number", "session", "medication", "run", "stim_id", "gvs_condition", "gvs_label", "block_position", "trial_in_run", "run_trial_index", "reward_level", "reward_high", "go_trial", "catch_trial", "excluded_session"]
    for measure in MEASURES:
        valid_col = f"valid_{measure.name.replace('+', '_plus_')}"
        piece = trials.loc[:, base_cols + [measure.output_column, valid_col]].copy()
        piece = piece.rename(columns={measure.output_column: "value", valid_col: "valid_measure"})
        piece["measure"] = measure.name
        piece["source_label"] = measure.source_label
        piece["unit"] = measure.unit
        piece["axis_label"] = measure.axis_label
        pieces.append(piece)
    return pd.concat(pieces, ignore_index=True)


def _successive_pairs(group, value_col):
    group = group.sort_values("trial_in_run")
    values = group[value_col].to_numpy(dtype=float)
    trial_idx = group["trial_in_run"].to_numpy(dtype=int)
    previous = values[:-1]
    following = values[1:]
    adjacent = np.diff(trial_idx) == 1
    keep = adjacent & np.isfinite(previous) & np.isfinite(following)
    return following[keep] - previous[keep]


def _rmssd_from_diffs(diffs):
    diffs = np.asarray(diffs, dtype=float)
    diffs = diffs[np.isfinite(diffs)]
    if diffs.size == 0:
        return np.nan, np.nan, 0
    mssd = float(np.mean(diffs**2))
    return float(np.sqrt(mssd)), mssd, int(diffs.size)


def _lag_autocorr_from_block(group, value_col, lag):
    group = group.sort_values("trial_in_run")
    values = group[value_col].to_numpy(dtype=float)
    trial_idx = group["trial_in_run"].to_numpy(dtype=int)
    if values.size <= lag:
        return np.nan
    x = values[:-lag]
    y = values[lag:]
    adjacent = (trial_idx[lag:] - trial_idx[:-lag]) == lag
    keep = adjacent & np.isfinite(x) & np.isfinite(y)
    if np.count_nonzero(keep) < 3:
        return np.nan
    return float(np.corrcoef(x[keep], y[keep])[0, 1])


def build_rt_block_metrics(trials):
    rows = []
    group_cols = ["subject", "subject_label", "session", "medication", "run", "stim_id", "gvs_condition", "gvs_label", "block_position", "excluded_session"]
    for key, group in trials.groupby(group_cols, sort=True):
        record = dict(zip(group_cols, key))
        go = group.loc[group["go_trial"]].copy()
        valid = go.loc[go["valid_rt_trial"]].copy()
        diffs = _successive_pairs(go, "RT_ms")
        rmssd, mssd, n_pairs = _rmssd_from_diffs(diffs)
        values = valid["RT_ms"].to_numpy(dtype=float)
        record.update(
            {
                "n_total_trials": int(len(group)),
                "n_go_trials": int(go.shape[0]),
                "n_valid_rt": int(values.size),
                "mean_rt_ms": float(np.nanmean(values)) if values.size else np.nan,
                "median_rt_ms": float(np.nanmedian(values)) if values.size else np.nan,
                "sd_rt_ms": float(np.nanstd(values, ddof=1)) if values.size > 1 else np.nan,
                "cv_rt": float(np.nanstd(values, ddof=1) / np.nanmean(values)) if values.size > 1 else np.nan,
                "rt_rmssd_ms": rmssd,
                "rt_mssd_ms2": mssd,
                "n_adjacent_pairs": n_pairs,
                "reward_high_fraction": float(valid["reward_high"].mean()) if values.size else np.nan,
                "lag1_autocorr": _lag_autocorr_from_block(go, "RT_ms", 1),
                "drift_slope_ms_per_trial": np.nan,
            }
        )
        if values.size >= 3:
            x = valid["trial_in_run"].to_numpy(dtype=float)
            y = valid["RT_ms"].to_numpy(dtype=float)
            if np.nanmax(x) > np.nanmin(x):
                record["drift_slope_ms_per_trial"] = float(np.polyfit(x, y, deg=1)[0])
        rows.append(record)
    return pd.DataFrame(rows).sort_values(group_cols).reset_index(drop=True)


def _measure_block_rmssd(trials, measure):
    rows = []
    group_cols = ["subject", "session", "medication", "run", "stim_id", "block_position", "excluded_session"]
    for key, group in trials.groupby(group_cols, sort=True):
        diffs = _successive_pairs(group.loc[group["go_trial"]], measure.output_column)
        rmssd, mssd, n_pairs = _rmssd_from_diffs(diffs)
        rows.append({**dict(zip(group_cols, key)), "measure": measure.name, "rmssd": rmssd, "mssd": mssd, "n_pairs": n_pairs, "unit": measure.unit})
    return pd.DataFrame(rows)


def build_feature_summary(trials, measure_long):
    rows = []
    included = measure_long.loc[~measure_long["excluded_session"] & measure_long["go_trial"]].copy()
    for measure in MEASURES:
        df = included.loc[included["measure"].eq(measure.name)].copy()
        valid = df.loc[df["valid_measure"] & np.isfinite(df["value"])].copy()
        missing_rate = 1.0 - (len(valid) / len(df)) if len(df) else np.nan

        run_summary = (valid.groupby(["subject", "session", "medication", "stim_id", "run"], as_index=False)["value"].mean())
        run_pivot = run_summary.pivot_table(index=["subject", "session", "medication", "stim_id"], columns="run", values="value")
        rel = _safe_pearson(run_pivot.get(1, pd.Series(dtype=float)), run_pivot.get(2, pd.Series(dtype=float)))

        med_summary = valid.groupby(["subject", "medication"], as_index=False)["value"].mean()
        med_pivot = med_summary.pivot(index="subject", columns="medication", values="value")
        med_delta = med_pivot["ON"] - med_pivot["OFF"] if {"OFF", "ON"}.issubset(med_pivot.columns) else pd.Series(dtype=float)
        med_stats = _safe_ttest_1samp(med_delta.to_numpy(dtype=float))

        reward_summary = (valid.loc[valid["reward_level"].isin(["low", "high"])] .groupby(["subject", "reward_level"], as_index=False)["value"] .mean())
        reward_pivot = reward_summary.pivot(index="subject", columns="reward_level", values="value")
        reward_delta = (reward_pivot["high"] - reward_pivot["low"] if {"low", "high"}.issubset(reward_pivot.columns) else pd.Series(dtype=float))
        reward_stats = _safe_ttest_1samp(reward_delta.to_numpy(dtype=float))

        stability = _measure_block_rmssd(trials, measure)
        stability = stability.loc[~stability["excluded_session"] & np.isfinite(stability["rmssd"])]
        stability_value = float(stability["rmssd"].mean()) if not stability.empty else np.nan

        if measure.name == "RT":
            reason = ("Kept: direct response-initiation latency, interpretable millisecond scale, " "available for most go trials, and compatible with neural trial timing.")
        elif measure.name == "PT":
            reason = "Excluded as primary: reflects pre-threshold pressure timing rather than response initiation."
        elif measure.name == "MT":
            reason = "Excluded as primary: movement execution duration is more biomechanical than vigour onset."
        elif measure.name == "RT+MT":
            reason = "Excluded as primary: composite mixes initiation and execution components."
        elif measure.name == "Vmax":
            reason = "Excluded as primary: force/velocity scale is less directly aligned with neural RT timing."
        else:
            reason = "Excluded as primary: peak force is a useful manipulation check but not the prespecified RT vigour measure."

        rows.append(
            {
                "Measure": measure.name,
                "Missing rate": missing_rate,
                "Test-retest reliability": f"run1-run2 r={_format_num(rel['r'], 2)} (n={rel['n']})",
                "Medication sensitivity": (f"ON-OFF={_format_num(med_stats['mean'], 2)} {measure.unit}; " f"{_format_p(float(med_stats['p']))}; n={med_stats['n']}"),
                "Reward sensitivity": (f"high-low={_format_num(reward_stats['mean'], 2)} {measure.unit}; " f"{_format_p(float(reward_stats['p']))}; n={reward_stats['n']}"),
                "Trial-wise stability": f"mean within-block RMSSD={_format_num(stability_value, 2)} {measure.unit}",
                "Reason kept/excluded": reason,
            }
        )
    return pd.DataFrame(rows)


def _fit_model(formula, data, group_col, model_name):
    if smf is None:
        return pd.DataFrame([{"model": model_name, "status": "skipped", "reason": "statsmodels unavailable"}])
    model_data = data.replace([np.inf, -np.inf], np.nan).dropna(subset=[formula.split("~")[0].strip(), group_col])
    rows = []
    fit = None
    fit_kind = "mixedlm_random_intercept_subject"
    reason = ""
    if model_data[group_col].nunique() >= 2:
        for method in ("powell", "bfgs", "cg", "nm", "lbfgs"):
            try:
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore")
                    fit = smf.mixedlm(formula, data=model_data, groups=model_data[group_col], re_formula="1").fit(reml=False, method=method, maxiter=300, disp=False)
                cov_re = float(fit.cov_re.iloc[0, 0]) if getattr(fit, "cov_re", pd.DataFrame()).size else np.nan
                finite_core = np.isfinite(getattr(fit, "llf", np.nan)) and np.all(np.isfinite(fit.params.to_numpy(dtype=float)))
                finite_key_bse = np.isfinite(float(fit.bse.iloc[0])) if len(fit.bse) else False
                if not finite_core or not finite_key_bse or not np.isfinite(cov_re) or cov_re <= 1e-10:
                    reason = f"singular or non-finite mixed model fit with {method}"
                    fit = None
                    continue
                break
            except Exception as exc:
                reason = str(exc)
                fit = None
    if fit is None:
        try:
            fit = smf.ols(formula, data=model_data).fit(cov_type="cluster", cov_kwds={"groups": model_data[group_col]})
            fit_kind = "ols_cluster_robust_by_subject_fallback"
        except Exception as exc:
            return pd.DataFrame([{"model": model_name, "status": "failed", "reason": reason or str(exc), "formula": formula}])

    conf = fit.conf_int()
    for term, coef in fit.params.items():
        rows.append(
            {
                "model": model_name,
                "status": "fit",
                "fit_kind": fit_kind,
                "term": term,
                "coefficient": float(coef),
                "std_error": float(fit.bse[term]) if term in fit.bse else np.nan,
                "z_or_t": float(fit.tvalues[term]) if term in fit.tvalues else np.nan,
                "p_value": float(fit.pvalues[term]) if term in fit.pvalues else np.nan,
                "ci95_low": float(conf.loc[term, 0]) if term in conf.index else np.nan,
                "ci95_high": float(conf.loc[term, 1]) if term in conf.index else np.nan,
                "n_observations": int(fit.nobs),
                "n_subjects": int(model_data[group_col].nunique()),
                "formula": formula,
            }
        )
    return pd.DataFrame(rows)


def build_model_tables(trials, block_metrics):
    rt_df = trials.loc[~trials["excluded_session"] & trials["valid_rt_trial"]].copy()
    rt_df = rt_df.loc[rt_df["reward_level"].isin(["low", "high"])].copy()
    rt_df["medication_on"] = (rt_df["medication"] == "ON").astype(int)
    rt_df["run"] = rt_df["run"].astype("category")
    rt_formula = ("RT_ms ~ medication_on * reward_high + medication_on * C(gvs_condition) " "+ C(run) + block_position + trial_in_run")

    var_df = block_metrics.loc[~block_metrics["excluded_session"] & np.isfinite(block_metrics["rt_rmssd_ms"])].copy()
    var_df["medication_on"] = (var_df["medication"] == "ON").astype(int)
    var_df["run"] = var_df["run"].astype("category")
    var_formula = "rt_rmssd_ms ~ medication_on + C(run) + block_position + reward_high_fraction + C(gvs_condition)"

    reward_rows = []
    for measure in [MEASURES[1], MEASURES[4], MEASURES[5]]:
        df = trials.loc[~trials["excluded_session"] & trials["go_trial"]].copy()
        valid_col = f"valid_{measure.name.replace('+', '_plus_')}"
        df = df.loc[df[valid_col] & df["reward_level"].isin(["low", "high"])].copy()
        df["value"] = df[measure.output_column]
        df["medication_on"] = (df["medication"] == "ON").astype(int)
        formula = "value ~ reward_high * medication_on + C(gvs_condition) + C(run) + block_position + trial_in_run"
        table = _fit_model(formula, df, "subject", f"reward_effect_{measure.name}")
        table["measure"] = measure.name
        table["unit"] = measure.unit
        reward_rows.append(table)

    return {"mean_rt_mixed_model": _fit_model(rt_formula, rt_df, "subject", "mean_rt_mixed_model"), "rt_variability_mixed_model": _fit_model(var_formula, var_df, "subject", "rt_variability_mixed_model"), "reward_mixed_models": pd.concat(reward_rows, ignore_index=True)}


def build_medication_deltas(trials, block_metrics):
    rt = (trials.loc[~trials["excluded_session"] & trials["valid_rt_trial"]] .groupby(["subject", "medication"], as_index=False)["RT_ms"] .mean())
    rt_pivot = rt.pivot(index="subject", columns="medication", values="RT_ms")
    var = (block_metrics.loc[~block_metrics["excluded_session"]] .groupby(["subject", "medication"], as_index=False)["rt_rmssd_ms"] .mean())
    var_pivot = var.pivot(index="subject", columns="medication", values="rt_rmssd_ms")
    subjects = sorted(set(rt_pivot.index) | set(var_pivot.index), key=_subject_number)
    rows = []
    for subject in subjects:
        off_rt = rt_pivot.loc[subject, "OFF"] if subject in rt_pivot.index and "OFF" in rt_pivot else np.nan
        on_rt = rt_pivot.loc[subject, "ON"] if subject in rt_pivot.index and "ON" in rt_pivot else np.nan
        off_var = var_pivot.loc[subject, "OFF"] if subject in var_pivot.index and "OFF" in var_pivot else np.nan
        on_var = var_pivot.loc[subject, "ON"] if subject in var_pivot.index and "ON" in var_pivot else np.nan
        rows.append({"subject": subject, "mean_rt_off_ms": off_rt, "mean_rt_on_ms": on_rt, "delta_rt_on_minus_off_ms": on_rt - off_rt if np.isfinite(off_rt) and np.isfinite(on_rt) else np.nan, "rt_rmssd_off_ms": off_var, "rt_rmssd_on_ms": on_var, "delta_rt_rmssd_on_minus_off_ms": on_var - off_var if np.isfinite(off_var) and np.isfinite(on_var) else np.nan})
    return pd.DataFrame(rows)


def build_reward_tables(trials, block_metrics):
    summary_rows = []
    for measure in [MEASURES[1], MEASURES[4], MEASURES[5]]:
        valid_col = f"valid_{measure.name.replace('+', '_plus_')}"
        df = trials.loc[~trials["excluded_session"] & trials["go_trial"] & trials[valid_col] & trials["reward_level"].isin(["low", "high"])].copy()
        grouped = (df.groupby(["subject", "medication", "reward_level"], as_index=False)[measure.output_column] .mean() .rename(columns={measure.output_column: "value"}))
        pivot = grouped.pivot_table(index=["subject", "medication"], columns="reward_level", values="value")
        pivot = pivot.reset_index()
        if {"low", "high"}.issubset(pivot.columns):
            pivot["high_minus_low"] = pivot["high"] - pivot["low"]
        pivot["measure"] = measure.name
        pivot["unit"] = measure.unit
        summary_rows.append(pivot)

    reward_var = block_metrics.loc[~block_metrics["excluded_session"] & np.isfinite(block_metrics["rt_rmssd_ms"])].copy()
    reward_var_summary = (reward_var.assign(reward_bin=pd.cut(reward_var["reward_high_fraction"], bins=[-0.01, 0.34, 0.67, 1.01])) .groupby(["subject", "medication", "reward_bin"], observed=False, as_index=False)["rt_rmssd_ms"] .mean())
    return {"reward_subject_session_summary": pd.concat(summary_rows, ignore_index=True), "reward_variability_by_high_reward_fraction": reward_var_summary}


def build_gvs_tables(trials, block_metrics, gvs_params):
    session_condition = (trials.loc[~trials["excluded_session"] & trials["valid_rt_trial"]] .groupby(["subject", "session", "medication", "stim_id", "gvs_condition", "gvs_label"], as_index=False)["RT_ms"] .mean() .rename(columns={"RT_ms": "mean_rt_ms"}))
    sham = session_condition.loc[session_condition["stim_id"].eq(1), ["subject", "session", "medication", "mean_rt_ms"]]
    sham = sham.rename(columns={"mean_rt_ms": "sham_rt_ms"})
    deltas = session_condition.loc[~session_condition["stim_id"].eq(1)].merge(sham, on=["subject", "session", "medication"], how="inner")
    deltas["delta_rt_ms_active_minus_sham"] = deltas["mean_rt_ms"] - deltas["sham_rt_ms"]

    rows = []
    params = gvs_params.set_index("stim_id")
    for stim_id, group in deltas.groupby("stim_id", sort=True):
        subject_avg = group.groupby("subject", as_index=False)["delta_rt_ms_active_minus_sham"].mean()
        test = _safe_ttest_1samp(subject_avg["delta_rt_ms_active_minus_sham"].to_numpy(dtype=float))
        off_mean = group.loc[group["medication"].eq("OFF"), "delta_rt_ms_active_minus_sham"].mean()
        on_mean = group.loc[group["medication"].eq("ON"), "delta_rt_ms_active_minus_sham"].mean()
        param = params.loc[stim_id] if stim_id in params.index else pd.Series(dtype=object)
        rows.append({"GVS condition": STIM_SHORT.get(int(stim_id), f"GVS{stim_id}"), "stim_id": int(stim_id), "Carrier": param.get("carrier_hz", np.nan), "Envelope": param.get("env_hz", np.nan), "family": param.get("family", ""), "main_freq_hz": param.get("main_freq_hz", np.nan), "RT change vs sham OFF": off_mean, "RT change vs sham ON": on_mean, "p-value": test["p"], "n_subjects": test["n"]})
    gvs_table = pd.DataFrame(rows)
    gvs_table["qFDR"] = _bh_fdr(gvs_table["p-value"].to_numpy(dtype=float))

    var_condition = (block_metrics.loc[~block_metrics["excluded_session"]] .groupby(["subject", "session", "medication", "stim_id", "gvs_condition", "gvs_label"], as_index=False)["rt_rmssd_ms"] .mean())
    var_sham = var_condition.loc[var_condition["stim_id"].eq(1), ["subject", "session", "medication", "rt_rmssd_ms"]]
    var_sham = var_sham.rename(columns={"rt_rmssd_ms": "sham_rt_rmssd_ms"})
    var_deltas = var_condition.loc[~var_condition["stim_id"].eq(1)].merge(var_sham, on=["subject", "session", "medication"], how="inner")
    var_deltas["delta_rt_rmssd_ms_active_minus_sham"] = var_deltas["rt_rmssd_ms"] - var_deltas["sham_rt_rmssd_ms"]
    return {"gvs_rt_session_deltas": deltas, "gvs_rt_variability_session_deltas": var_deltas, "gvs_effect_summary_table": gvs_table}


def build_temporal_tables(trials, block_metrics):
    included = block_metrics.loc[~block_metrics["excluded_session"]].copy()
    rows = []
    for label, group in [("overall", included)] + [(med, g) for med, g in included.groupby("medication", sort=True)]:
        rows.append(
            {
                "group": label,
                "n_blocks": int(group.shape[0]),
                "mean_rt_sd_ms": float(group["sd_rt_ms"].mean()),
                "mean_rt_cv": float(group["cv_rt"].mean()),
                "mean_rt_rmssd_ms": float(group["rt_rmssd_ms"].mean()),
                "mean_rt_mssd_ms2": float(group["rt_mssd_ms2"].mean()),
                "mean_lag1_autocorrelation": float(group["lag1_autocorr"].mean()),
                "mean_drift_slope_ms_per_trial": float(group["drift_slope_ms_per_trial"].mean()),
                "n_adjacent_pairs": int(group["n_adjacent_pairs"].sum()),
                "variability_boundary_rule": "successive differences computed within subject-session-run-GVS blocks only",
            }
        )
    temporal_summary = pd.DataFrame(rows)

    acf_rows = []
    for lag in range(1, 6):
        values = []
        for _, group in trials.loc[~trials["excluded_session"] & trials["go_trial"]].groupby(["subject", "session", "run", "stim_id"], sort=False):
            values.append(_lag_autocorr_from_block(group, "RT_ms", lag))
        values = np.asarray(values, dtype=float)
        acf_rows.append({"lag": lag, "mean_autocorrelation": float(np.nanmean(values)), "sem_autocorrelation": float(stats.sem(values[np.isfinite(values)])) if np.count_nonzero(np.isfinite(values)) > 1 else np.nan, "n_blocks": int(np.count_nonzero(np.isfinite(values)))})
    acf_table = pd.DataFrame(acf_rows)
    return {"temporal_structure_summary": temporal_summary, "rt_acf_by_lag": acf_table}


def build_neural_residual_sensitivity(trials, neural_mean_path, neural_variability_path):
    valid = trials.loc[~trials["excluded_session"] & trials["valid_rt_trial"]].copy()
    valid = valid.loc[valid["reward_level"].isin(["low", "high"])].copy()
    if smf is None or valid.empty:
        return {"reward_residualized_rt_trials": pd.DataFrame(), "neural_residual_sensitivity": pd.DataFrame()}

    formula = "RT_ms ~ reward_high + C(medication) + C(gvs_condition) + C(run) + block_position"
    fit = smf.ols(formula, data=valid).fit()
    valid["RT_residual_ms"] = fit.resid
    residual_trials = valid.loc[:, ["subject", "subject_label", "session", "medication", "run", "stim_id", "gvs_label", "block_position", "trial_in_run", "RT_ms", "RT_residual_ms", "reward_high"]].copy()

    residual_block_rows = []
    for key, group in residual_trials.groupby(["subject_label", "subject", "session", "medication", "run", "stim_id", "gvs_label"], sort=True):
        diffs = _successive_pairs(group, "RT_residual_ms")
        rmssd, mssd, n_pairs = _rmssd_from_diffs(diffs)
        residual_block_rows.append({"subject_label": key[0], "subject": key[1], "session": key[2], "medication": key[3], "run": key[4], "stim_id": key[5], "gvs_label": key[6], "mean_rt_residual_ms": float(group["RT_residual_ms"].mean()), "rt_residual_rmssd_ms": rmssd, "rt_residual_mssd_ms2": mssd, "n_residual_pairs": n_pairs})
    residual_block = pd.DataFrame(residual_block_rows)
    residual_session = (residual_block.groupby(["subject_label", "subject", "session", "medication", "stim_id", "gvs_label"], as_index=False) .agg(mean_rt_residual_ms=("mean_rt_residual_ms", "mean"), rt_residual_rmssd_ms=("rt_residual_rmssd_ms", "mean"), n_residual_pairs=("n_residual_pairs", "sum")))
    sham = residual_session.loc[residual_session["stim_id"].eq(1), ["subject_label", "session", "medication", "mean_rt_residual_ms", "rt_residual_rmssd_ms"]].rename(columns={"mean_rt_residual_ms": "sham_mean_rt_residual_ms", "rt_residual_rmssd_ms": "sham_rt_residual_rmssd_ms"})
    residual_delta = residual_session.loc[~residual_session["stim_id"].eq(1)].merge(sham, on=["subject_label", "session", "medication"], how="inner")
    residual_delta["delta_mean_rt_residual_ms"] = (residual_delta["mean_rt_residual_ms"] - residual_delta["sham_mean_rt_residual_ms"])
    residual_delta["delta_rt_residual_rmssd_ms"] = (residual_delta["rt_residual_rmssd_ms"] - residual_delta["sham_rt_residual_rmssd_ms"])

    sensitivity_rows = []
    if neural_mean_path.exists():
        neural_mean = pd.read_csv(neural_mean_path)
        neural_session = (neural_mean.groupby(["subject", "session", "medication", "stim_id", "stim_short"], as_index=False)["mean_proj"] .mean())
        neural_sham = neural_session.loc[neural_session["stim_id"].eq(1), ["subject", "session", "medication", "mean_proj"]].rename(columns={"mean_proj": "sham_mean_proj"})
        neural_delta = neural_session.loc[~neural_session["stim_id"].eq(1)].merge(neural_sham, on=["subject", "session", "medication"], how="inner")
        neural_delta["delta_mean_proj"] = neural_delta["mean_proj"] - neural_delta["sham_mean_proj"]
        merged = neural_delta.merge(residual_delta, left_on=["subject", "session", "medication", "stim_id"], right_on=["subject_label", "session", "medication", "stim_id"], how="inner")
        corr = _safe_pearson(merged["delta_mean_proj"], merged["delta_mean_rt_residual_ms"])
        sensitivity_rows.append({"analysis": "mean_projection_delta_vs_reward_residualized_mean_rt_delta", "n": corr["n"], "r": corr["r"], "p_value": corr["p"], "interpretation": "tests whether sham-referenced neural projection mean is explained by residual RT"})
    if neural_variability_path.exists():
        neural_var = pd.read_csv(neural_variability_path)
        merged = neural_var.loc[~neural_var["stim_id"].eq(1)].merge(residual_delta, left_on=["subject", "session", "medication", "stim_id"], right_on=["subject_label", "session", "medication", "stim_id"], how="inner")
        corr = _safe_pearson(merged["delta_consecutive_diff2"], merged["delta_rt_residual_rmssd_ms"])
        sensitivity_rows.append({"analysis": "projection_variability_delta_vs_reward_residualized_rt_rmssd_delta", "n": corr["n"], "r": corr["r"], "p_value": corr["p"], "interpretation": "tests whether sham-referenced neural projection variability tracks residual RT variability"})
    return {"reward_residualized_rt_trials": residual_trials, "reward_residualized_rt_session_deltas": residual_delta, "neural_residual_sensitivity": pd.DataFrame(sensitivity_rows)}


def build_subject_summary(trials, block_metrics, reward_tables, gvs_tables, hand_table):
    med = build_medication_deltas(trials, block_metrics)
    reward = reward_tables["reward_subject_session_summary"]
    reward_rt = reward.loc[reward["measure"].eq("RT")].groupby("subject", as_index=False)["high_minus_low"].mean()
    reward_rt = reward_rt.rename(columns={"high_minus_low": "reward_effect_high_minus_low_rt_ms"})
    gvs = gvs_tables["gvs_rt_session_deltas"].groupby("subject", as_index=False)["delta_rt_ms_active_minus_sham"].mean()
    gvs = gvs.rename(columns={"delta_rt_ms_active_minus_sham": "gvs_sensitivity_active_minus_sham_rt_ms"})
    summary = med.merge(reward_rt, on="subject", how="outer").merge(gvs, on="subject", how="outer")
    summary = summary.merge(hand_table, on="subject", how="left")
    flow_status = (trials.groupby(["subject", "session"], as_index=False)["excluded_session"] .any() .query("excluded_session") .groupby("subject")["session"] .apply(lambda x: ", ".join(str(int(v)) for v in sorted(x))) .reset_index(name="excluded_sessions"))
    summary = summary.merge(flow_status, on="subject", how="left")
    summary["notes"] = ""
    summary.loc[summary["excluded_sessions"].notna(), "notes"] = ("session " + summary.loc[summary["excluded_sessions"].notna(), "excluded_sessions"].astype(str) + " excluded")
    summary.loc[summary["non_dominant_task_hand"].fillna(False), "notes"] = summary.loc[summary["non_dominant_task_hand"].fillna(False), "notes"].map(lambda x: _combine_notes([x, "non-dominant task hand"]))
    summary["_sort"] = summary["subject"].map(_subject_number)
    summary = summary.sort_values("_sort").drop(columns="_sort").reset_index(drop=True)
    return summary


def build_distribution_summary(trials):
    valid = trials.loc[~trials["excluded_session"] & trials["valid_rt_trial"]].copy()
    rows = []
    groups = [("overall", []), ("medication", ["medication"]), ("run", ["run"]), ("block", ["block_position"]), ("trial_position", ["trial_in_run"]), ("subject", ["subject"])]
    for label, cols in groups:
        if not cols:
            iterable = [(("overall",), valid)]
            cols = ["level"]
        else:
            iterable = valid.groupby(cols, sort=True)
        for key, group in iterable:
            key_tuple = key if isinstance(key, tuple) else (key,)
            record = {"summary_type": label}
            for col, value in zip(cols, key_tuple):
                record[col] = value
            values = group["RT_ms"].to_numpy(dtype=float)
            record.update({"n": int(values.size), "mean_rt_ms": float(np.mean(values)), "median_rt_ms": float(np.median(values)), "sd_rt_ms": float(np.std(values, ddof=1)) if values.size > 1 else np.nan, "q25_rt_ms": float(np.quantile(values, 0.25)), "q75_rt_ms": float(np.quantile(values, 0.75))})
            rows.append(record)
    return pd.DataFrame(rows)


def plot_measure_distributions(measure_long, out_dir):
    valid = measure_long.loc[~measure_long["excluded_session"] & measure_long["valid_measure"]].copy()
    with plt.rc_context(PAPER_STYLE):
        fig, axes = plt.subplots(2, 3, figsize=(11.0, 6.4))
        axes = axes.ravel()
        for ax, measure in zip(axes, MEASURES):
            df = valid.loc[valid["measure"].eq(measure.name)]
            values = df["value"].to_numpy(dtype=float)
            values = values[np.isfinite(values)]
            ax.hist(values, bins=35, color="#4C78A8", alpha=0.78, density=True)
            ax.set_title(measure.name)
            ax.set_xlabel(measure.axis_label)
            ax.set_ylabel("Density")
        fig.suptitle("Behavioural feature distributions", y=1.01)
        fig.tight_layout()
        _save_figure(fig, out_dir / "measure_distribution_all_features")


def plot_rt_medication_distribution(trials, out_dir):
    valid = trials.loc[~trials["excluded_session"] & trials["valid_rt_trial"]].copy()
    colors = {"OFF": "#4C78A8", "ON": "#D55E00"}
    with plt.rc_context(BOLD_PAPER_STYLE):
        fig, ax = plt.subplots(figsize=(7.4, 4.8))
        x_grid = np.linspace(valid["RT_ms"].quantile(0.005), valid["RT_ms"].quantile(0.995), 400)
        for med, group in valid.groupby("medication", sort=True):
            values = group["RT_ms"].to_numpy(dtype=float)
            values = values[np.isfinite(values)]
            if values.size < 2:
                continue
            density = stats.gaussian_kde(values)(x_grid)
            color = colors.get(med, "0.3")
            ax.fill_between(x_grid, density, color=color, alpha=0.26, label=med)
            ax.plot(x_grid, density, color=color, linewidth=2.0)
        ax.set_xlabel("RT (ms)")
        ax.set_ylabel("Density")
        ax.grid(False)
        ax.legend(frameon=False)
        _bold_figure_text(fig)
        fig.tight_layout()
        _save_figure(fig, out_dir / "supp_fig_behaviour_1_rt_distribution_medication")


def plot_rt_paired_medication(subject_summary, out_dir):
    with plt.rc_context(BOLD_PAPER_STYLE):
        fig, ax = plt.subplots(figsize=(4.25, 5.15))
        x = np.array([0.0, 0.55])
        paired = subject_summary.dropna(subset=["mean_rt_off_ms", "mean_rt_on_ms"])
        for _, row in paired.iterrows():
            line_color = "#111827" if row.get("non_dominant_task_hand", False) else "#8CA0B3"
            dot_color = "#111827" if row.get("non_dominant_task_hand", False) else "#5F7488"
            lw = 1.8 if row.get("non_dominant_task_hand", False) else 1.0
            ax.plot(x, [row["mean_rt_off_ms"], row["mean_rt_on_ms"]], color=line_color, alpha=0.55, linewidth=lw)
            ax.scatter(x, [row["mean_rt_off_ms"], row["mean_rt_on_ms"]], color=dot_color, s=24, zorder=3)
        missing = subject_summary.loc[subject_summary[["mean_rt_off_ms", "mean_rt_on_ms"]].isna().any(axis=1)]
        for _, row in missing.iterrows():
            if np.isfinite(row.get("mean_rt_off_ms", np.nan)):
                ax.scatter([x[0]], [row["mean_rt_off_ms"]], marker="x", color="#C23B4B", s=56, linewidth=2.0)
            if np.isfinite(row.get("mean_rt_on_ms", np.nan)):
                ax.scatter([x[1]], [row["mean_rt_on_ms"]], marker="x", color="#C23B4B", s=56, linewidth=2.0)
        means = paired[["mean_rt_off_ms", "mean_rt_on_ms"]].mean().to_numpy(dtype=float)
        sems = paired[["mean_rt_off_ms", "mean_rt_on_ms"]].sem().to_numpy(dtype=float)
        ax.errorbar(x, means, yerr=sems, color="#1F2937", marker="D", markersize=6, linewidth=2.2, capsize=4)
        ax.set_xticks(x)
        ax.set_xticklabels(["OFF", "ON"])
        ax.set_xlim(-0.10, 0.65)
        ax.set_ylabel("Mean RT (ms)")
        _bold_figure_text(fig)
        fig.tight_layout()
        _save_figure(fig, out_dir / "supp_fig_behaviour_2_subject_paired_rt_off_on")


def plot_rt_run_block_trial(trials, out_dir):
    valid = trials.loc[~trials["excluded_session"] & trials["valid_rt_trial"]].copy()
    with plt.rc_context(BOLD_PAPER_STYLE):
        fig, axes = plt.subplots(1, 3, figsize=(14.6, 4.35), gridspec_kw={"width_ratios": [0.65, 1.2, 1.35]})
        panels = [("run", "Run", [1, 2], "#0072B2"), ("block_position", "Block position in run", list(range(1, 10)), "#009E73"), ("trial_in_run", "Trial position within block", list(range(1, 11)), "#D55E00")]
        for ax, (col, xlabel, order, color) in zip(axes, panels):
            summary = valid.groupby(col)["RT_ms"].agg(["mean", "count", "std"]).reindex(order).reset_index()
            sem = summary["std"] / np.sqrt(summary["count"])
            ax.errorbar(summary[col], summary["mean"], yerr=sem, marker="o", color=color, capsize=3)
            ax.set_xlabel(xlabel)
            ax.set_ylabel("RT (ms)")
            ax.set_xticks(order)
        axes[0].set_xlim(0.5, 2.5)
        axes[0].set_title("Run")
        axes[1].set_title("Block")
        axes[2].set_title("Trial")
        _bold_figure_text(fig)
        fig.tight_layout(pad=0.35, w_pad=2.2)
        _save_figure(fig, out_dir / "supp_fig_behaviour_3_rt_by_run_block_trial")


def plot_trial_flow_qc_stacked_bar(trial_flow, out_dir):
    df = trial_flow.copy()
    df["subject_number"] = df["Subject"].astype(str).str.extract(r"(\d+)").astype(int)
    subject_order = df.sort_values("subject_number")[["Subject", "subject_number"]].drop_duplicates("Subject")
    subjects = subject_order["Subject"].tolist()
    subject_labels = [f"P{int(subject_number):03d}" for subject_number in subject_order["subject_number"]]
    medications = [med for med in ["OFF", "ON"] if med in set(df["Medication"])]
    colors = {"Valid RT": "#4C78A8", "Missed": "#FF7F0E", "Invalid RT": "#B07AA1", "Excluded": "#9CA3AF"}
    x = np.arange(len(subjects), dtype=float)
    width = 0.31
    offsets = {1: -0.18, 2: 0.18}

    with plt.rc_context(BOLD_PAPER_STYLE):
        fig, axes = plt.subplots(len(medications), 1, figsize=(15.2, 8.0), sharex=True, sharey=True)
        axes = np.atleast_1d(axes)
        for ax, medication in zip(axes, medications, strict=True):
            med_df = df.loc[df["Medication"].eq(medication)]
            for run, offset in offsets.items():
                run_df = med_df.loc[med_df["Run"].eq(run)].set_index("Subject").reindex(subjects)
                xpos = x + offset
                included = run_df["Session status"].fillna("").eq("included").to_numpy()
                valid = run_df["Valid RT trials"].fillna(0).to_numpy(dtype=float)
                missed = run_df["Missed trials"].fillna(0).to_numpy(dtype=float)
                invalid = run_df["Invalid RT"].fillna(0).to_numpy(dtype=float)
                ax.bar(xpos[included], valid[included], width=width, color=colors["Valid RT"], edgecolor="white", linewidth=0.5)
                ax.bar(xpos[included], missed[included], width=width, bottom=valid[included], color=colors["Missed"], edgecolor="white", linewidth=0.5)
                ax.bar(xpos[included], invalid[included], width=width, bottom=valid[included] + missed[included], color=colors["Invalid RT"], edgecolor="white", linewidth=0.5)
                excluded = ~included & run_df["Go trials"].notna().to_numpy()
                ax.bar(xpos[excluded], run_df["Go trials"].fillna(0).to_numpy(dtype=float)[excluded], width=width, color=colors["Excluded"], edgecolor="white", linewidth=0.5)
            ax.axhline(81, color="#111827", linestyle=(0, (4, 3)), linewidth=2.2)
            ax.set_ylim(0, 86)
            ax.set_yticks([0, 20, 40, 60, 80])
            ax.grid(axis="y", color="#E5E7EB", linewidth=0.8)
        axes[-1].set_xticks(x)
        axes[-1].set_xticklabels(subject_labels, rotation=45, ha="right")
        fig.supxlabel("Subject", fontsize=20, fontweight="bold", y=0.02)
        fig.supylabel("Number of trials", fontsize=20, fontweight="bold", x=0.02)
        legend_handles = [Patch(facecolor=colors["Valid RT"], label="Valid RT"), Patch(facecolor=colors["Missed"], label="Missed"), Patch(facecolor=colors["Invalid RT"], label="Invalid RT")]
        fig.legend(handles=legend_handles, loc="upper center", ncol=3, frameon=False, bbox_to_anchor=(0.5, 1.0))
        _bold_figure_text(fig)
        fig.tight_layout(rect=(0.04, 0.04, 1.0, 0.94), h_pad=1.6)
        _save_figure(fig, out_dir / "supp_fig_behaviour_4_trial_flow_qc_stacked_bar")


def plot_medication_delta(subject_summary, out_dir):
    with plt.rc_context(PAPER_STYLE):
        fig, axes = plt.subplots(1, 2, figsize=(8.8, 4.2))
        specs = [("delta_rt_on_minus_off_ms", "Delta RT=ON - OFF (ms)"), ("delta_rt_rmssd_on_minus_off_ms", "Delta RT RMSSD=ON - OFF (ms)")]
        for ax, (col, label) in zip(axes, specs):
            df = subject_summary.loc[np.isfinite(subject_summary[col])].copy()
            x = np.arange(len(df))
            colors = np.where(df["non_dominant_task_hand"].fillna(False), "#111827", "#4C78A8")
            ax.axhline(0, color="0.35", linestyle="--", linewidth=0.9)
            ax.scatter(x, df[col], color=colors, s=28)
            mean, ci_low, ci_high = _mean_ci(df[col])
            ax.errorbar([len(df) + 0.7], [mean], yerr=[[mean - ci_low], [ci_high - mean]], fmt="D", color="#D55E00")
            ax.set_xticks(x)
            ax.set_xticklabels(df["subject"].str.replace("PSPD", "P", regex=False), rotation=90)
            ax.set_ylabel(label)
        fig.tight_layout()
        _save_figure(fig, out_dir / "medication_delta_rt_and_variability_subjects")


def plot_reward(reward_tables, block_metrics, out_dir):
    reward = reward_tables["reward_subject_session_summary"]
    with plt.rc_context(PAPER_STYLE):
        fig, axes = plt.subplots(2, 2, figsize=(10.4, 7.0))
        for ax, measure in zip(axes.flat[:3], ["RT", "Vmax", "Pmax"]):
            df = reward.loc[reward["measure"].eq(measure)].copy()
            long = df.melt(id_vars=["subject", "medication", "measure", "unit"], value_vars=["low", "high"], var_name="reward_level", value_name="value")
            x_base = {"low": 0.0, "high": 1.0}
            for med, offset, color in [("OFF", -0.08, "#4C78A8"), ("ON", 0.08, "#D55E00")]:
                med_df = long.loc[long["medication"].eq(med)]
                means = med_df.groupby("reward_level")["value"].mean()
                sems = med_df.groupby("reward_level")["value"].sem()
                xs = np.array([x_base["low"] + offset, x_base["high"] + offset])
                ax.errorbar(xs, [means.get("low", np.nan), means.get("high", np.nan)], yerr=[sems.get("low", np.nan), sems.get("high", np.nan)], marker="o", color=color, capsize=3, label=med)
            ax.set_xticks([0, 1])
            ax.set_xticklabels(["Low", "High"])
            ax.set_title(measure)
            unit = df["unit"].dropna().iloc[0] if not df.empty else ""
            ax.set_ylabel(f"{measure} ({unit})")
        var = block_metrics.loc[~block_metrics["excluded_session"]].copy()
        ax = axes.flat[3]
        ax.scatter(var["reward_high_fraction"], var["rt_rmssd_ms"], s=18, alpha=0.35, color="#009E73")
        if var["reward_high_fraction"].nunique() > 1:
            x = var["reward_high_fraction"].to_numpy(dtype=float)
            y = var["rt_rmssd_ms"].to_numpy(dtype=float)
            keep = np.isfinite(x) & np.isfinite(y)
            slope, intercept = np.polyfit(x[keep], y[keep], deg=1)
            x_fit = np.linspace(0, 1, 50)
            ax.plot(x_fit, slope * x_fit + intercept, color="#111827", linewidth=1.5)
        ax.set_xlabel("High-reward fraction in block")
        ax.set_ylabel("RT RMSSD (ms)")
        ax.set_title("RT variability")
        axes.flat[0].legend(frameon=False)
        fig.suptitle("Reward manipulation checks", y=1.01)
        fig.tight_layout()
        _save_figure(fig, out_dir / "reward_manipulation_check")


def plot_gvs(gvs_tables, block_metrics, out_dir):
    deltas = gvs_tables["gvs_rt_session_deltas"]
    var_deltas = gvs_tables["gvs_rt_variability_session_deltas"]
    with plt.rc_context(PAPER_STYLE):
        fig, axes = plt.subplots(2, 2, figsize=(11.0, 7.2))
        rt_summary = (block_metrics.loc[~block_metrics["excluded_session"]] .groupby(["stim_id", "gvs_label"], as_index=False)["mean_rt_ms"] .agg(["mean", "sem"]) .reset_index())
        axes[0, 0].errorbar(rt_summary["stim_id"], rt_summary["mean"], yerr=rt_summary["sem"], marker="o", color="#4C78A8")
        axes[0, 0].set_xticks(rt_summary["stim_id"])
        axes[0, 0].set_xticklabels(rt_summary["gvs_label"], rotation=35, ha="right")
        axes[0, 0].set_ylabel("Mean RT (ms)")
        axes[0, 0].set_title("RT by GVS")

        var_summary = (block_metrics.loc[~block_metrics["excluded_session"]] .groupby(["stim_id", "gvs_label"], as_index=False)["rt_rmssd_ms"] .agg(["mean", "sem"]) .reset_index())
        axes[0, 1].errorbar(var_summary["stim_id"], var_summary["mean"], yerr=var_summary["sem"], marker="o", color="#009E73")
        axes[0, 1].set_xticks(var_summary["stim_id"])
        axes[0, 1].set_xticklabels(var_summary["gvs_label"], rotation=35, ha="right")
        axes[0, 1].set_ylabel("RT RMSSD (ms)")
        axes[0, 1].set_title("RT variability by GVS")

        rt_delta_summary = deltas.groupby(["stim_id", "gvs_label"], as_index=False)["delta_rt_ms_active_minus_sham"].agg(["mean", "sem"]).reset_index()
        axes[1, 0].axhline(0, color="0.35", linestyle="--", linewidth=0.9)
        axes[1, 0].errorbar(rt_delta_summary["stim_id"], rt_delta_summary["mean"], yerr=rt_delta_summary["sem"], marker="o", color="#D55E00")
        axes[1, 0].set_xticks(rt_delta_summary["stim_id"])
        axes[1, 0].set_xticklabels(rt_delta_summary["gvs_label"], rotation=35, ha="right")
        axes[1, 0].set_ylabel("Active - sham RT (ms)")
        axes[1, 0].set_title("Sham-referenced RT")

        for med, color in [("OFF", "#4C78A8"), ("ON", "#D55E00")]:
            med_df = deltas.loc[deltas["medication"].eq(med)]
            summary = med_df.groupby(["stim_id", "gvs_label"], as_index=False)["delta_rt_ms_active_minus_sham"].agg(["mean", "sem"]).reset_index()
            axes[1, 1].errorbar(summary["stim_id"], summary["mean"], yerr=summary["sem"], marker="o", color=color, label=med)
        axes[1, 1].axhline(0, color="0.35", linestyle="--", linewidth=0.9)
        axes[1, 1].set_xticks(rt_delta_summary["stim_id"])
        axes[1, 1].set_xticklabels(rt_delta_summary["gvs_label"], rotation=35, ha="right")
        axes[1, 1].set_ylabel("Active - sham RT (ms)")
        axes[1, 1].set_title("By medication")
        axes[1, 1].legend(frameon=False)
        fig.suptitle("GVS behavioural analysis", y=1.01)
        fig.tight_layout()
        _save_figure(fig, out_dir / "gvs_behaviour_summary")

        fig, ax = plt.subplots(figsize=(6.5, 4.0))
        var_summary = var_deltas.groupby(["stim_id", "gvs_label"], as_index=False)["delta_rt_rmssd_ms_active_minus_sham"].agg(["mean", "sem"]).reset_index()
        ax.axhline(0, color="0.35", linestyle="--", linewidth=0.9)
        ax.errorbar(var_summary["stim_id"], var_summary["mean"], yerr=var_summary["sem"], marker="o", color="#009E73")
        ax.set_xticks(var_summary["stim_id"])
        ax.set_xticklabels(var_summary["gvs_label"], rotation=35, ha="right")
        ax.set_ylabel("Active - sham RT RMSSD (ms)")
        ax.set_title("Sham-referenced RT variability by GVS")
        fig.tight_layout()
        _save_figure(fig, out_dir / "gvs_rt_variability_sham_referenced")


def plot_temporal(temporal_tables, trials, out_dir):
    acf = temporal_tables["rt_acf_by_lag"]
    valid = trials.loc[~trials["excluded_session"] & trials["valid_rt_trial"]].copy()
    drift = valid.groupby("run_trial_index")["RT_ms"].agg(["mean", "sem"]).reset_index()
    with plt.rc_context(PAPER_STYLE):
        fig, axes = plt.subplots(1, 2, figsize=(9.0, 3.8))
        axes[0].errorbar(acf["lag"], acf["mean_autocorrelation"], yerr=acf["sem_autocorrelation"], marker="o", color="#4C78A8")
        axes[0].axhline(0, color="0.35", linestyle="--", linewidth=0.9)
        axes[0].set_xlabel("Lag")
        axes[0].set_ylabel("Autocorrelation")
        axes[0].set_title("RT autocorrelation")
        axes[1].errorbar(drift["run_trial_index"], drift["mean"], yerr=drift["sem"], color="#D55E00", linewidth=1.5)
        axes[1].set_xlabel("Task trial index within run")
        axes[1].set_ylabel("RT (ms)")
        axes[1].set_title("RT drift")
        fig.tight_layout()
        _save_figure(fig, out_dir / "temporal_structure_rt_acf_and_drift")


def plot_subject_summary(subject_summary, out_dir):
    df = subject_summary.copy()
    df["subject_short"] = df["subject"].str.replace("PSPD", "P", regex=False)
    y = np.arange(len(df))
    with plt.rc_context(PAPER_STYLE):
        fig, axes = plt.subplots(1, 4, figsize=(13.5, max(5.5, 0.32 * len(df) + 1.8)), sharey=True)
        axes[0].scatter(df["mean_rt_off_ms"], y, color="#4C78A8", label="OFF", s=22)
        axes[0].scatter(df["mean_rt_on_ms"], y, color="#D55E00", label="ON", s=22)
        for yi, row in zip(y, df.itertuples(index=False)):
            if np.isfinite(row.mean_rt_off_ms) and np.isfinite(row.mean_rt_on_ms):
                axes[0].plot([row.mean_rt_off_ms, row.mean_rt_on_ms], [yi, yi], color="0.70", linewidth=0.8, zorder=0)
        axes[0].set_xlabel("Mean RT (ms)")
        axes[0].legend(frameon=False, fontsize=8)

        axes[1].scatter(df["rt_rmssd_off_ms"], y, color="#4C78A8", s=22)
        axes[1].scatter(df["rt_rmssd_on_ms"], y, color="#D55E00", s=22)
        for yi, row in zip(y, df.itertuples(index=False)):
            if np.isfinite(row.rt_rmssd_off_ms) and np.isfinite(row.rt_rmssd_on_ms):
                axes[1].plot([row.rt_rmssd_off_ms, row.rt_rmssd_on_ms], [yi, yi], color="0.70", linewidth=0.8, zorder=0)
        axes[1].set_xlabel("RT RMSSD (ms)")

        axes[2].axvline(0, color="0.35", linestyle="--", linewidth=0.8)
        axes[2].scatter(df["reward_effect_high_minus_low_rt_ms"], y, color="#009E73", s=24)
        axes[2].set_xlabel("High - low reward RT\n(ms; negative = faster)")

        axes[3].axvline(0, color="0.35", linestyle="--", linewidth=0.8)
        axes[3].scatter(df["gvs_sensitivity_active_minus_sham_rt_ms"], y, color="#7C3AED", s=24)
        axes[3].set_xlabel("GVS sensitivity\nactive-sham RT (ms)")

        for ax in axes:
            ax.set_ylim(-0.8, len(df) - 0.2)
            ax.invert_yaxis()
        axes[0].set_yticks(y)
        labels = []
        for _, row in df.iterrows():
            label = row["subject_short"]
            notes = []
            if bool(row.get("non_dominant_task_hand", False)):
                notes.append("ND")
            if isinstance(row.get("excluded_sessions", np.nan), str):
                notes.append("S" + row["excluded_sessions"] + " excl.")
            labels.append(label + (" (" + ", ".join(notes) + ")" if notes else ""))
        axes[0].set_yticklabels(labels, fontsize=8)
        for ax in axes[1:]:
            ax.tick_params(axis="y", labelleft=False)
        fig.suptitle("Subject-level behavioural heterogeneity", y=1.005)
        fig.tight_layout()
        _save_figure(fig, out_dir / "subject_level_behaviour_summary")


def write_definitions(out_dir, args, unresolved_hand):
    text = [
        "# Behaviour Supplement Assumptions and Definitions",
        "",
        "## Source Files",
        "",
        f"- Behaviour metrics: `{args.metrics_dir}`",
        f"- Consolidated metadata: `{args.consolidated_dir}`",
        f"- GVS order inventory: `{args.gvs_order}`",
        f"- GVS parameter metadata: `{args.gvs_params}`",
        "",
        "The Windows source path supplied by the user maps to `/mnt/TeamShare/Data_Masterfile/.../PRECISIONSTIM_PD_Data_Results/Behaviour/Behaviour_metrics_revised` in this environment.",
        "",
        "## Trial Definitions",
        "",
        "- Total trials: all 20 entries per GVS condition split into run 1 trials 1-10 and run 2 trials 11-20.",
        "- Catch trials: trials with `res.catchtrials == 1`; this was cross-checked against reward code -5.",
        "- Go trials: all non-catch trials.",
        "- Missed trials: go trials marked as no squeeze in `res.sqrwd_nosq == 1`.",
        "- Invalid RT trials: go trials that were not missed but did not have finite positive source `1/RT`.",
        "- Valid RT trials: go trials with finite positive source `1/RT`; RT is reported as `1000 / (1/RT)` in ms.",
        "",
        "Very-late responses with finite positive RT were retained as valid RT trials to match the existing behavioural scripts, which use finite positive RT rather than reward-success filtering.",
        "",
        "## Exclusions",
        "",
        "- PSPD017 session 1 is marked and excluded from group analyses because one recording run was excessively noisy.",
        "- The primary analysis uses the session numbering in the GVS/fMRI inventory: OFF=session 1 and ON=session 2.",
        "- The older `figures/reward_effects` output used a subject-specific PSPD017 medication/session override, so reward effects are recomputed here from raw files for consistency.",
        "",
        "## RT Variability",
        "",
        "RT variability is root mean squared successive difference (RMSSD) in milliseconds. It is computed only for adjacent valid RT pairs within the same subject, session, run, and GVS block. Transitions across GVS blocks, rest periods, and run boundaries are excluded.",
    ]
    if unresolved_hand:
        text.extend(
            [
                "",
                "## Non-Dominant Task Hand",
                "",
                "The prompt notes that one right-handed participant used the left hand because dominant-hand tremor was severe. The local behavioural metrics, consolidated metadata, GVS order files, and available demographics tables did not identify which subject this was. The script therefore exposes `--non-dominant-task-hand-subject PSPD###` so the exact participant can be marked when the ID is known.",
            ]
        )
    (out_dir / "assumptions_and_definitions.md").write_text("\n".join(text) + "\n", encoding="utf-8")


def write_methods_results(out_dir, trials, trial_flow, feature_summary, med_deltas, gvs_summary, temporal_summary, neural_sensitivity):
    included_valid = trials.loc[~trials["excluded_session"] & trials["valid_rt_trial"]]
    n_subjects = included_valid["subject"].nunique()
    n_sessions = included_valid.groupby(["subject", "session"]).ngroups
    n_valid = included_valid.shape[0]
    rt_mean = included_valid["RT_ms"].mean()
    rt_sd = included_valid["RT_ms"].std()
    med_rt = _safe_ttest_1samp(med_deltas["delta_rt_on_minus_off_ms"].to_numpy(dtype=float))
    med_var = _safe_ttest_1samp(med_deltas["delta_rt_rmssd_on_minus_off_ms"].to_numpy(dtype=float))
    rt_feature = feature_summary.loc[feature_summary["Measure"].eq("RT")].iloc[0]
    gvs_min = gvs_summary.sort_values("qFDR", na_position="last").iloc[0] if not gvs_summary.empty else None
    temporal = temporal_summary.loc[temporal_summary["group"].eq("overall")].iloc[0]

    lines = [
        "# Behavioural Supplementary Methods and Results",
        "",
        "## Methods",
        "",
        ("Behavioural analyses used the first six columns of each revised metrics file: " "1/PT, 1/RT, 1/MT, 1/(RT+MT), Vmax, and Pmax. The first four columns were " "reciprocal time measures and were converted to milliseconds for reporting. " "Vmax and Pmax were kept on their raw feature scales."),
        "",
        ("Each GVS condition contributed 20 trials per session. Trials 1-10 were assigned " "to run 1 and trials 11-20 to run 2. Actual block position within each run was " "read from `data/gvs_order_by_subject_session_run.tsv`."),
        "",
        ("Catch trials were identified from `res.catchtrials` and verified against reward code -5. " "Go trials were non-catch trials. Missed trials were go trials marked as no-squeeze " "in `res.sqrwd_nosq`. Valid RT trials were go trials with finite positive source 1/RT; " "invalid RT trials were non-missed go trials without finite positive 1/RT."),
        "",
        ("RT variability was defined as RMSSD in milliseconds, computed from adjacent valid " "RT pairs within subject-session-run-GVS blocks. This excludes transitions across " "GVS blocks, rest gaps, and run boundaries."),
        "",
        "## Results",
        "",
        (f"After excluding PSPD017 session 1, the analysis retained {n_valid} valid RT trials " f"from {n_sessions} sessions and {n_subjects} subjects. Mean RT was {rt_mean:.1f} ms " f"(SD across trials {rt_sd:.1f} ms)."),
        "",
        (
            "RT was retained as the primary behavioural vigour measure because it is a direct, "
            "interpretable response-initiation latency with low missingness and a stable mapping "
            "to the trial-wise neural analyses. In the feature comparison table, RT missing rate "
            f"was {float(rt_feature['Missing rate']):.3f}; {rt_feature['Test-retest reliability']}. "
            "Force features were retained as manipulation checks rather than the primary vigour outcome."
        ),
        "",
        (f"Medication changed mean RT by ON-OFF={med_rt['mean']:.1f} ms " f"(95% CI {med_rt['ci95_low']:.1f}, {med_rt['ci95_high']:.1f}; " f"{_format_p(float(med_rt['p']))}; n={med_rt['n']}). " f"Medication changed RT RMSSD by ON-OFF={med_var['mean']:.1f} ms " f"(95% CI {med_var['ci95_low']:.1f}, {med_var['ci95_high']:.1f}; " f"{_format_p(float(med_var['p']))}; n={med_var['n']})."),
    ]
    if gvs_min is not None:
        lines.extend(["", (f"GVS effects on RT were small. The smallest FDR-adjusted condition-level " f"RT effect was {gvs_min['GVS condition']} with q={gvs_min['qFDR']:.3f} " f"and {_format_p(float(gvs_min['p-value']))}.")])
    lines.extend(["", (f"Across blocks, mean RT RMSSD was {temporal['mean_rt_rmssd_ms']:.1f} ms, " f"mean lag-1 autocorrelation was {temporal['mean_lag1_autocorrelation']:.3f}, " f"and mean linear drift was {temporal['mean_drift_slope_ms_per_trial']:.2f} ms/trial.")])
    if not neural_sensitivity.empty:
        lines.append("")
        for row in neural_sensitivity.itertuples(index=False):
            lines.append(f"Reward-residualized neural sensitivity ({row.analysis}) gave r={_format_num(row.r, 3)}, " f"{_format_p(float(row.p_value))}, n={int(row.n)}.")
    lines.extend(["", "The generated CSV tables contain the full trial flow, feature comparison, model coefficients, GVS condition tests, temporal structure, and subject-level heterogeneity summaries."])
    (out_dir / "supplementary_methods_results_behaviour.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def copy_reused_manifest(args, out_dir):
    rows = []
    for label, path in [("prior_gvs_behaviour", args.prior_gvs_dir), ("prior_reward_behaviour", args.prior_reward_dir), ("neural_mean_projection", args.neural_mean_table), ("neural_projection_variability", args.neural_variability_table)]:
        rows.append({"resource": label, "path": str(path), "exists": path.exists(), "use": "read for context/sensitivity"})
    manifest = pd.DataFrame(rows)
    _write_csv(manifest, out_dir / "reused_outputs_manifest.csv")
    for source in [args.prior_gvs_dir / "gvs_behaviour_report.md", args.prior_reward_dir / "reward_rt_report.md"]:
        if source.exists():
            shutil.copy2(source, out_dir / f"legacy_{source.name}")
    return manifest


def main():
    args = _parse_args()
    out_dir = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    non_dom_subject = _pspd_id(args.non_dominant_task_hand_subject) if args.non_dominant_task_hand_subject else None
    unresolved_hand = non_dom_subject is None
    gvs_params = _load_gvs_params(args.gvs_params)
    hand_table = _load_dominant_hand(args.dominant_hand_table, non_dom_subject)
    order_map = _load_order_map(args.gvs_order)
    copy_reused_manifest(args, out_dir)

    trials = build_trial_table(args.metrics_dir, args.consolidated_dir, gvs_params, order_map, hand_table)
    measure_long = build_measure_long(trials)
    trial_flow = build_trial_flow(trials)
    block_metrics = build_rt_block_metrics(trials)
    feature_summary = build_feature_summary(trials, measure_long)
    med_deltas = build_medication_deltas(trials, block_metrics)
    reward_tables = build_reward_tables(trials, block_metrics)
    gvs_tables = build_gvs_tables(trials, block_metrics, gvs_params)
    temporal_tables = build_temporal_tables(trials, block_metrics)
    neural_tables = build_neural_residual_sensitivity(trials, args.neural_mean_table, args.neural_variability_table)
    subject_summary = build_subject_summary(trials, block_metrics, reward_tables, gvs_tables, hand_table)
    distribution_summary = build_distribution_summary(trials)

    _write_csv(trials, out_dir / "behaviour_supp_trial_table.csv")
    _write_csv(measure_long, out_dir / "behaviour_supp_measure_long.csv")
    _write_csv(trial_flow, out_dir / "table_behaviour_trial_flow_qc.csv")
    _write_markdown_table(trial_flow, out_dir / "table_behaviour_trial_flow_qc.md")
    _write_csv(block_metrics, out_dir / "rt_block_variability_metrics.csv")
    _write_csv(feature_summary, out_dir / "table_behaviour_feature_comparison.csv")
    _write_markdown_table(feature_summary, out_dir / "table_behaviour_feature_comparison.md")
    _write_csv(med_deltas, out_dir / "medication_subject_deltas.csv")
    _write_csv(reward_tables["reward_subject_session_summary"], out_dir / "reward_subject_session_summary.csv")
    _write_csv(reward_tables["reward_variability_by_high_reward_fraction"], out_dir / "reward_variability_by_high_reward_fraction.csv")
    for name, df in gvs_tables.items():
        _write_csv(df, out_dir / f"{name}.csv")
    for name, df in temporal_tables.items():
        _write_csv(df, out_dir / f"{name}.csv")
    for name, df in neural_tables.items():
        _write_csv(df, out_dir / f"{name}.csv")
    _write_csv(subject_summary, out_dir / "subject_level_behaviour_summary.csv")
    _write_csv(distribution_summary, out_dir / "rt_distribution_summary.csv")
    _write_csv(gvs_params, out_dir / "gvs_condition_parameter_table.csv")

    if not args.skip_models:
        model_tables = build_model_tables(trials, block_metrics)
        for name, df in model_tables.items():
            _write_csv(df, out_dir / f"{name}.csv")
    else:
        model_tables = {}

    plot_measure_distributions(measure_long, out_dir)
    plot_rt_medication_distribution(trials, out_dir)
    plot_rt_paired_medication(subject_summary, out_dir)
    plot_rt_run_block_trial(trials, out_dir)
    plot_trial_flow_qc_stacked_bar(trial_flow, out_dir)
    plot_medication_delta(subject_summary, out_dir)
    plot_reward(reward_tables, block_metrics, out_dir)
    plot_gvs(gvs_tables, block_metrics, out_dir)
    plot_temporal(temporal_tables, trials, out_dir)
    plot_subject_summary(subject_summary, out_dir)

    write_definitions(out_dir, args, unresolved_hand)
    write_methods_results(out_dir, trials, trial_flow, feature_summary, med_deltas, gvs_tables["gvs_effect_summary_table"], temporal_tables["temporal_structure_summary"], neural_tables["neural_residual_sensitivity"])

    manifest = {
        "n_raw_trial_rows": int(trials.shape[0]),
        "n_included_valid_rt_trials": int((~trials["excluded_session"] & trials["valid_rt_trial"]).sum()),
        "n_subjects_included": int(trials.loc[~trials["excluded_session"], "subject"].nunique()),
        "excluded_subject_sessions": [{"subject": subject, "session": session, "reason": reason} for (subject, session), reason in EXCLUDED_SUBJECT_SESSIONS.items()],
        "non_dominant_task_hand_subject": non_dom_subject,
        "outputs_directory": str(out_dir),
    }
    (out_dir / "analysis_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
