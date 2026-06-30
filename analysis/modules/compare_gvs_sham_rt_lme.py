#!/usr/bin/env python3
"""Compare trial-level RT across GVS conditions with mixed-effects models.

By default this script reads revised behav_metrics .mat files, where cells are
ordered sham, GVS1, GVS2, ..., GVS8. It reads the selected RT metric column,
z-scores 1/RT within each subject/session, and tests condition effects
separately for each medication/session using trial-level LME models. RT in ms is
kept only as a derived descriptive column.
"""

from __future__ import annotations

import argparse
import json
import re
import warnings
from pathlib import Path
from typing import Any

warnings.filterwarnings("ignore", message="Unable to import Axes3D.*")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats
from scipy.io import loadmat
from statsmodels.formula.api import mixedlm
from statsmodels.stats.multitest import fdrcorrection


ROOT = Path(__file__).resolve().parents[2]
SCRIPT_DIR = Path(__file__).resolve().parent
DEFAULT_METRICS_DIR = Path(
    "/mnt/TeamShare/Data_Masterfile/H20-00572_All-Dressed/"
    "AllDressed_WorkOnData/Sepideh/Behaviour_metrics_revised"
)
FALLBACK_METRICS_DIR = Path(
    "/mnt/TeamShare/Data_Masterfile/H20-00572_All-Dressed/"
    "PRECISIONSTIM_PD_Data_Results/Behaviour/Behaviour_metrics_revised"
)
DEFAULT_BEHAVIOUR_DIR = ROOT / "data" / "external" / "behaviour_npy"
DEFAULT_INVENTORY = ROOT / "data" / "processed" / "gvs_connectivity" / "common" / "run_condition_inventory.csv"
DEFAULT_OUT_DIR = ROOT / "results" / "main" / "figure_08_gvs_effects" / "rt_lme"

CONDITION_ORDER = ["sham"] + [f"GVS{i}" for i in range(1, 9)]
ALPHA = 0.05
MIXEDLM_METHODS = ("powell", "nm", "lbfgs", "bfgs", "cg")
SESSION_BY_MEDICATION = {"OFF": 1, "ON": 2}
METRICS_RE = re.compile(r"^(PSPD\d+)_(OFF|ON)_behav_metrics\.mat$")

STYLE = {
    "font.family": "sans-serif",
    "font.sans-serif": ["Liberation Sans", "Arial", "DejaVu Sans"],
    "font.size": 11,
    "axes.labelsize": 12,
    "axes.titlesize": 12,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "axes.spines.top": False,
    "axes.spines.right": False,
    "axes.grid": True,
    "grid.color": "#E5E7EB",
    "grid.linewidth": 0.8,
    "grid.alpha": 0.8,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--data-source",
        choices=("metrics_mat", "npy"),
        default="metrics_mat",
        help="Use revised .mat metrics files, or the older local .npy RT vectors.",
    )
    parser.add_argument(
        "--metrics-dir",
        type=coerce_data_path,
        default=DEFAULT_METRICS_DIR,
        help="Directory containing PSPD###_OFF/ON_behav_metrics.mat files.",
    )
    parser.add_argument(
        "--rt-column",
        type=int,
        default=2,
        help="One-based RT metric column in each behav_metrics cell. Default 2.",
    )
    parser.add_argument("--behaviour-dir", type=Path, default=DEFAULT_BEHAVIOUR_DIR)
    parser.add_argument("--inventory", type=Path, default=DEFAULT_INVENTORY)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument(
        "--input-kind",
        choices=("inverse_rt", "rt_ms", "rt_s"),
        default="inverse_rt",
        help="How to interpret values in the .npy RT arrays.",
    )
    parser.add_argument(
        "--no-run-fixed-effect",
        action="store_true",
        help="Fit condition-only LME models instead of adjusting for run.",
    )
    parser.add_argument(
        "--use-inventory-trial-ranges",
        action="store_true",
        help=(
            "Use trial_start/trial_stop from the inventory. By default the script "
            "uses stim_order to map all 90 local RT trials into nine 10-trial blocks."
        ),
    )
    parser.add_argument("--alpha", type=float, default=ALPHA)
    return parser.parse_args()


def coerce_data_path(value: str | Path) -> Path:
    text = str(value)
    if re.match(r"^[A-Za-z]:\\", text):
        drive, rest = text[:2], text[3:]
        if drive.upper() == "Z:":
            return Path("/mnt/TeamShare") / rest.replace("\\", "/")
    return Path(text)


def pspd_id(subject: str) -> str:
    match = re.search(r"(\d+)", str(subject))
    if match is None:
        raise ValueError(f"Could not parse subject number from {subject!r}")
    return f"PSPD{int(match.group(1)):03d}"


def to_rt_ms(values: np.ndarray, input_kind: str) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    out = np.full(values.shape, np.nan, dtype=float)
    valid = np.isfinite(values) & (values > 0)
    if input_kind == "inverse_rt":
        out[valid] = 1000.0 / values[valid]
    elif input_kind == "rt_s":
        out[valid] = 1000.0 * values[valid]
    else:
        out[valid] = values[valid]
    return out


def to_inverse_rt(values: np.ndarray, input_kind: str) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    out = np.full(values.shape, np.nan, dtype=float)
    valid = np.isfinite(values) & (values > 0)
    if input_kind == "inverse_rt":
        out[valid] = values[valid]
    elif input_kind == "rt_s":
        out[valid] = 1.0 / values[valid]
    else:
        out[valid] = 1000.0 / values[valid]
    return out


def select_metric_column(values: np.ndarray, rt_column_one_based: int, path: Path) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    if values.ndim == 1:
        return values
    if values.ndim != 2:
        raise RuntimeError(f"{path} must be a 1D vector or 2D metric matrix, got shape {values.shape}")

    rt_column_index = int(rt_column_one_based) - 1
    if rt_column_index < 0:
        raise ValueError("--rt-column must be one-based and positive")
    if values.shape[1] <= rt_column_index:
        raise RuntimeError(
            f"{path} has {values.shape[1]} columns; cannot read one-based column {rt_column_one_based}"
        )
    return values[:, rt_column_index]


def add_inverse_rt_zscores(trials: pd.DataFrame) -> pd.DataFrame:
    trials = trials.copy()
    trials["inverse_rt_z"] = trials.groupby(["subject", "session"], group_keys=False)["inverse_rt"].transform(
        zscore_series
    )
    # Compatibility alias for the existing model/plot code.
    trials["rt_z"] = trials["inverse_rt_z"]
    return trials


def condition_from_stim(stim_id: int) -> tuple[str, str]:
    code = f"gvs-{int(stim_id):02d}"
    label = "sham" if int(stim_id) == 1 else f"GVS{int(stim_id) - 1}"
    return code, label


def parse_stim_order(value: Any) -> list[int]:
    return [int(part.strip()) for part in str(value).split(",") if str(part).strip()]


def parse_metrics_filename(path: Path) -> tuple[str, str, int]:
    match = METRICS_RE.match(path.name)
    if match is None:
        raise ValueError(f"Unrecognized metrics filename: {path.name}")
    subject = match.group(1)
    medication = match.group(2).upper()
    return subject, medication, SESSION_BY_MEDICATION[medication]


def subpd_id(subject: str) -> str:
    match = re.search(r"(\d+)", str(subject))
    if match is None:
        raise ValueError(f"Could not parse subject number from {subject!r}")
    return f"sub-pd{int(match.group(1)):03d}"


def load_metric_cells(path: Path) -> list[np.ndarray]:
    mat = loadmat(path, squeeze_me=False, struct_as_record=False)
    if "behav_metrics" not in mat:
        raise RuntimeError(f"{path} does not contain behav_metrics")
    cells = np.asarray(mat["behav_metrics"], dtype=object).ravel()
    metrics = [np.asarray(cell, dtype=float) for cell in cells]
    if len(metrics) < len(CONDITION_ORDER):
        raise RuntimeError(f"{path} has {len(metrics)} metric cells; expected at least {len(CONDITION_ORDER)}")
    return metrics[: len(CONDITION_ORDER)]


def load_trial_table_from_metrics(metrics_dir: Path, input_kind: str, rt_column_one_based: int) -> pd.DataFrame:
    if not metrics_dir.exists() and metrics_dir == DEFAULT_METRICS_DIR and FALLBACK_METRICS_DIR.exists():
        metrics_dir = FALLBACK_METRICS_DIR
    metric_paths = sorted(metrics_dir.glob("PSPD*_behav_metrics.mat"))
    if not metric_paths:
        raise FileNotFoundError(f"No PSPD*_behav_metrics.mat files found under {metrics_dir}")

    rt_column_index = int(rt_column_one_based) - 1
    if rt_column_index < 0:
        raise ValueError("--rt-column must be one-based and positive")

    rows: list[dict[str, Any]] = []
    for metric_path in metric_paths:
        subject, medication, session = parse_metrics_filename(metric_path)
        subject_label = subpd_id(subject)
        metric_cells = load_metric_cells(metric_path)
        for stim_index, metric in enumerate(metric_cells, start=1):
            if metric.ndim != 2:
                raise RuntimeError(f"{metric_path} cell {stim_index} has shape {metric.shape}; expected 2D")
            if metric.shape[1] <= rt_column_index:
                raise RuntimeError(
                    f"{metric_path} cell {stim_index} has {metric.shape[1]} columns; "
                    f"cannot read one-based column {rt_column_one_based}"
                )
            condition_code, condition_label = condition_from_stim(stim_index)
            source_values = metric[:, rt_column_index].astype(float)
            inverse_rt = to_inverse_rt(source_values, input_kind)
            rt_ms = to_rt_ms(source_values, input_kind)
            for trial_in_condition, (source_value, inverse_rt_value, rt_value) in enumerate(
                zip(source_values, inverse_rt, rt_ms, strict=True),
                start=1,
            ):
                run = 1 if trial_in_condition <= 10 else 2
                rows.append(
                    {
                        "subject": subject_label,
                        "source_subject": subject,
                        "session": session,
                        "medication": medication,
                        "run": run,
                        "condition_code": condition_code,
                        "condition_label": condition_label,
                        "active_gvs": condition_code != "gvs-01",
                        "trial_in_condition": int(trial_in_condition),
                        "task_trial_index_zero_based": int(trial_in_condition - 1),
                        "source_value": float(source_value) if np.isfinite(source_value) else np.nan,
                        "inverse_rt": float(inverse_rt_value) if np.isfinite(inverse_rt_value) else np.nan,
                        "rt_ms": float(rt_value) if np.isfinite(rt_value) else np.nan,
                        "behaviour_path": str(metric_path),
                    }
                )

    trials = pd.DataFrame(rows)
    trials["condition_label"] = pd.Categorical(trials["condition_label"], categories=CONDITION_ORDER, ordered=True)
    return add_inverse_rt_zscores(trials)


def load_trial_table(
    behaviour_dir: Path,
    inventory_path: Path,
    input_kind: str,
    rt_column_one_based: int,
) -> pd.DataFrame:
    inventory = pd.read_csv(inventory_path)
    rows: list[dict[str, Any]] = []
    cache: dict[Path, np.ndarray] = {}

    for _, block in inventory.iterrows():
        subject = str(block["subject"])
        session = int(block["session"])
        run = int(block["run"])
        path = behaviour_dir / f"{pspd_id(subject)}_ses_{session}_run_{run}.npy"
        if path not in cache:
            if not path.exists():
                warnings.warn(f"Missing behaviour file: {path}", RuntimeWarning)
                continue
            cache[path] = select_metric_column(np.load(path, allow_pickle=False), rt_column_one_based, path)

        values = cache[path]
        start = int(block["trial_start"])
        stop = int(block["trial_stop"])
        if start < 0 or stop > len(values) or start >= stop:
            warnings.warn(
                f"Skipping invalid trial range {start}:{stop} for {path.name}",
                RuntimeWarning,
            )
            continue

        source_values = values[start:stop]
        inverse_rt = to_inverse_rt(source_values, input_kind)
        rt_ms = to_rt_ms(source_values, input_kind)
        for offset, (source_value, inverse_rt_value, rt_value) in enumerate(
            zip(source_values, inverse_rt, rt_ms, strict=True),
            start=1,
        ):
            rows.append(
                {
                    "subject": subject,
                    "session": session,
                    "medication": str(block["medication"]).upper(),
                    "run": run,
                    "condition_code": str(block["condition_code"]),
                    "condition_label": str(block["condition_label"]),
                    "active_gvs": str(block["condition_code"]) != "gvs-01",
                    "trial_in_condition": offset,
                    "task_trial_index_zero_based": start + offset - 1,
                    "source_value": float(source_value) if np.isfinite(source_value) else np.nan,
                    "inverse_rt": float(inverse_rt_value) if np.isfinite(inverse_rt_value) else np.nan,
                    "rt_ms": float(rt_value) if np.isfinite(rt_value) else np.nan,
                    "behaviour_path": str(path),
                }
            )

    if not rows:
        raise RuntimeError("No trial rows were loaded. Check --behaviour-dir and --inventory.")

    trials = pd.DataFrame(rows)
    trials["condition_label"] = pd.Categorical(trials["condition_label"], categories=CONDITION_ORDER, ordered=True)
    return add_inverse_rt_zscores(trials)


def load_trial_table_from_stim_order(
    behaviour_dir: Path,
    inventory_path: Path,
    input_kind: str,
    rt_column_one_based: int,
) -> pd.DataFrame:
    inventory = pd.read_csv(inventory_path)
    required = {"subject", "session", "medication", "run", "stim_order"}
    missing = sorted(required - set(inventory.columns))
    if missing:
        raise RuntimeError(f"{inventory_path} is missing required columns: {', '.join(missing)}")

    rows: list[dict[str, Any]] = []
    cache: dict[Path, np.ndarray] = {}
    runs = inventory.drop_duplicates(["subject", "session", "medication", "run"]).copy()
    runs = runs.sort_values(["subject", "session", "run"]).reset_index(drop=True)

    for run_row in runs.itertuples(index=False):
        subject = str(run_row.subject)
        session = int(run_row.session)
        run = int(run_row.run)
        path = behaviour_dir / f"{pspd_id(subject)}_ses_{session}_run_{run}.npy"
        if path not in cache:
            if not path.exists():
                warnings.warn(f"Missing behaviour file: {path}", RuntimeWarning)
                continue
            cache[path] = select_metric_column(np.load(path, allow_pickle=False), rt_column_one_based, path)

        values = cache[path]
        stim_order = parse_stim_order(run_row.stim_order)
        if len(stim_order) == 0:
            warnings.warn(f"Skipping empty stim_order for {subject} session {session} run {run}", RuntimeWarning)
            continue

        trials_per_block = len(values) // len(stim_order)
        if trials_per_block <= 0:
            warnings.warn(f"Skipping too-short behaviour file for {path.name}", RuntimeWarning)
            continue

        for block_position, stim_id in enumerate(stim_order, start=1):
            start = (block_position - 1) * trials_per_block
            stop = start + trials_per_block
            condition_code, condition_label = condition_from_stim(stim_id)
            block_values = values[start:stop]
            inverse_rt = to_inverse_rt(block_values, input_kind)
            rt_ms = to_rt_ms(block_values, input_kind)
            for offset, (source_value, inverse_rt_value, rt_value) in enumerate(
                zip(block_values, inverse_rt, rt_ms, strict=True),
                start=1,
            ):
                rows.append(
                    {
                        "subject": subject,
                        "session": session,
                        "medication": str(run_row.medication).upper(),
                        "run": run,
                        "condition_code": condition_code,
                        "condition_label": condition_label,
                        "active_gvs": condition_code != "gvs-01",
                        "block_position": int(block_position),
                        "trial_in_condition": int(offset),
                        "task_trial_index_zero_based": int(start + offset - 1),
                        "source_value": float(source_value) if np.isfinite(source_value) else np.nan,
                        "inverse_rt": float(inverse_rt_value) if np.isfinite(inverse_rt_value) else np.nan,
                        "rt_ms": float(rt_value) if np.isfinite(rt_value) else np.nan,
                        "behaviour_path": str(path),
                    }
                )

    if not rows:
        raise RuntimeError("No trial rows were loaded. Check --behaviour-dir and --inventory.")

    trials = pd.DataFrame(rows)
    trials["condition_label"] = pd.Categorical(trials["condition_label"], categories=CONDITION_ORDER, ordered=True)
    return add_inverse_rt_zscores(trials)


def zscore_series(values: pd.Series) -> pd.Series:
    mean = values.mean(skipna=True)
    sd = values.std(skipna=True, ddof=0)
    if not np.isfinite(sd) or sd == 0:
        return pd.Series(np.nan, index=values.index)
    return (values - mean) / sd


def fit_mixed_model(formula: str, data: pd.DataFrame) -> tuple[Any, str, list[str]]:
    last_error: Exception | None = None
    for method in MIXEDLM_METHODS:
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            try:
                model = mixedlm(formula, data=data, groups=data["subject"], re_formula="1")
                result = model.fit(reml=False, method=method, maxiter=600, disp=False)
            except Exception as exc:  # pragma: no cover - depends on optimizer state
                last_error = exc
                continue
        warning_text = [str(item.message) for item in caught]
        if np.isfinite(result.llf):
            return result, method, warning_text
    raise RuntimeError(f"MixedLM failed for formula {formula!r}: {last_error}")


def condition_term(result: Any, condition_label: str) -> str | None:
    suffix = f"[T.{condition_label}]"
    matches = [term for term in result.params.index if str(term).startswith("C(condition_label") and str(term).endswith(suffix)]
    return matches[0] if matches else None


def fit_session_models(
    trials: pd.DataFrame,
    *,
    include_run: bool,
    alpha: float,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    anova_rows: list[dict[str, Any]] = []
    pairwise_rows: list[dict[str, Any]] = []
    term_rows: list[dict[str, Any]] = []
    model_trials = trials.loc[np.isfinite(trials["inverse_rt_z"])].copy()

    run_term = " + C(run)" if include_run else ""
    full_formula = f'inverse_rt_z ~ C(condition_label, Treatment(reference="sham")){run_term}'
    null_formula = "inverse_rt_z ~ C(run)" if include_run else "inverse_rt_z ~ 1"

    group_order = (
        model_trials[["session", "medication"]]
        .drop_duplicates()
        .sort_values(["session", "medication"])
        .itertuples(index=False, name=None)
    )
    for session, medication in group_order:
        data = model_trials.loc[
            model_trials["session"].eq(session) & model_trials["medication"].eq(medication)
        ].copy()
        full, full_method, full_warnings = fit_mixed_model(full_formula, data)
        null, null_method, null_warnings = fit_mixed_model(null_formula, data)

        lr_stat = max(0.0, 2.0 * (float(full.llf) - float(null.llf)))
        df_diff = int(full.df_modelwc - null.df_modelwc)
        omnibus_p = float(stats.chi2.sf(lr_stat, df_diff)) if df_diff > 0 else np.nan
        common = {
            "session": int(session),
            "medication": medication,
            "n_trials": int(len(data)),
            "n_subjects": int(data["subject"].nunique()),
        }
        anova_rows.append(
            {
                **common,
                "test": "mixedlm_likelihood_ratio_anova",
                "full_formula": full_formula,
                "null_formula": null_formula,
                "chi2": lr_stat,
                "df": df_diff,
                "p_value": omnibus_p,
                "full_llf": float(full.llf),
                "null_llf": float(null.llf),
                "full_aic": float(full.aic),
                "null_aic": float(null.aic),
                "full_optimizer": full_method,
                "null_optimizer": null_method,
                "full_converged": bool(getattr(full, "converged", False)),
                "null_converged": bool(getattr(null, "converged", False)),
                "full_warnings": " | ".join(dict.fromkeys(full_warnings)),
                "null_warnings": " | ".join(dict.fromkeys(null_warnings)),
            }
        )

        for term in full.params.index:
            term_rows.append(
                {
                    **common,
                    "term": str(term),
                    "coef": float(full.params.get(term, np.nan)),
                    "se": float(full.bse.get(term, np.nan)),
                    "z_stat": float(full.tvalues.get(term, np.nan)),
                    "p_value": float(full.pvalues.get(term, np.nan)),
                }
            )

        p_values = []
        new_rows = []
        for condition in CONDITION_ORDER[1:]:
            term = condition_term(full, condition)
            condition_data = data.loc[data["condition_label"].eq(condition)]
            p_value = float(full.pvalues.get(term, np.nan)) if term is not None else np.nan
            row = {
                **common,
                "condition_label": condition,
                "reference": "sham",
                "term": term,
                "coef_z_inverse_rt": float(full.params.get(term, np.nan)) if term is not None else np.nan,
                "se": float(full.bse.get(term, np.nan)) if term is not None else np.nan,
                "z_stat": float(full.tvalues.get(term, np.nan)) if term is not None else np.nan,
                "p_value": p_value,
                "direction": (
                    "faster_than_sham"
                    if np.isfinite(p_value) and float(full.params.get(term, np.nan)) > 0
                    else "slower_than_sham"
                ),
                "n_trials_condition": int(len(condition_data)),
                "mean_inverse_rt_z_condition": float(condition_data["inverse_rt_z"].mean()),
                "median_inverse_rt_z_condition": float(condition_data["inverse_rt_z"].median()),
            }
            new_rows.append(row)
            p_values.append(p_value)

        finite = np.isfinite(p_values)
        q_values = np.full(len(p_values), np.nan, dtype=float)
        rejected = np.full(len(p_values), False, dtype=bool)
        if np.any(finite):
            rejected_finite, q_finite = fdrcorrection(np.asarray(p_values, dtype=float)[finite], alpha=alpha)
            q_values[finite] = q_finite
            rejected[finite] = rejected_finite
        for row, q_value, is_rejected in zip(new_rows, q_values, rejected, strict=True):
            row["q_fdr_within_session_medication"] = float(q_value) if np.isfinite(q_value) else np.nan
            row["significant_fdr"] = bool(is_rejected)
            row["significant_raw_p"] = bool(np.isfinite(row["p_value"]) and row["p_value"] < alpha)
            pairwise_rows.append(row)

    return pd.DataFrame(anova_rows), pd.DataFrame(pairwise_rows), pd.DataFrame(term_rows)


def summarize_conditions(trials: pd.DataFrame) -> pd.DataFrame:
    finite = trials.loc[np.isfinite(trials["inverse_rt_z"])].copy()
    trial_summary = (
        finite.groupby(["session", "medication", "condition_code", "condition_label"], observed=True)
        .agg(
            n_trials=("inverse_rt_z", "count"),
            n_subjects=("subject", "nunique"),
            mean_inverse_rt=("inverse_rt", "mean"),
            median_inverse_rt=("inverse_rt", "median"),
            mean_rt_ms=("rt_ms", "mean"),
            median_rt_ms=("rt_ms", "median"),
            mean_inverse_rt_z=("inverse_rt_z", "mean"),
            median_inverse_rt_z=("inverse_rt_z", "median"),
            sem_inverse_rt_z_trial=("inverse_rt_z", "sem"),
        )
        .reset_index()
    )
    subject_condition = (
        finite.groupby(["session", "medication", "subject", "condition_code", "condition_label"], observed=True)
        .agg(
            subject_mean_inverse_rt=("inverse_rt", "mean"),
            subject_mean_rt_ms=("rt_ms", "mean"),
            subject_mean_inverse_rt_z=("inverse_rt_z", "mean"),
        )
        .reset_index()
    )
    subject_summary = (
        subject_condition.groupby(["session", "medication", "condition_code", "condition_label"], observed=True)
        .agg(
            mean_subject_inverse_rt=("subject_mean_inverse_rt", "mean"),
            sem_subject_inverse_rt=("subject_mean_inverse_rt", "sem"),
            mean_subject_rt_ms=("subject_mean_rt_ms", "mean"),
            sem_subject_rt_ms=("subject_mean_rt_ms", "sem"),
            mean_subject_inverse_rt_z=("subject_mean_inverse_rt_z", "mean"),
            sem_subject_inverse_rt_z=("subject_mean_inverse_rt_z", "sem"),
        )
        .reset_index()
    )
    summary = trial_summary.merge(
        subject_summary,
        on=["session", "medication", "condition_code", "condition_label"],
        how="left",
    )
    summary["_condition_order"] = summary["condition_label"].map({name: i for i, name in enumerate(CONDITION_ORDER)})
    summary = summary.sort_values(["session", "medication", "_condition_order"]).drop(columns="_condition_order")
    return summary


def format_p(value: float) -> str:
    if not np.isfinite(value):
        return "p = NA"
    if value < 0.001:
        return "p < 0.001"
    return f"p = {value:.3f}"


def condition_display_label(condition: str, condition_code: str | None = None) -> str:
    if str(condition).lower() == "sham" or str(condition_code).lower() == "gvs-01":
        return "Sham"
    if re.fullmatch(r"GVS\d+", str(condition)):
        return str(condition)
    code_match = re.search(r"gvs-(\d+)", str(condition_code).lower()) if condition_code is not None else None
    if code_match:
        return f"GVS{int(code_match.group(1)) - 1}"
    return str(condition)


def condition_code_order(trials: pd.DataFrame) -> list[str]:
    codes = sorted(
        trials["condition_code"].dropna().astype(str).unique(),
        key=lambda value: int(re.search(r"(\d+)$", value).group(1)) if re.search(r"(\d+)$", value) else 999,
    )
    return codes


def paired_t_pvalue(active: np.ndarray, reference: np.ndarray) -> float:
    active = np.asarray(active, dtype=float)
    reference = np.asarray(reference, dtype=float)
    keep = np.isfinite(active) & np.isfinite(reference)
    if np.count_nonzero(keep) < 2:
        return np.nan
    diff = active[keep] - reference[keep]
    if np.nanstd(diff, ddof=1) == 0:
        return 1.0 if np.nanmean(diff) == 0 else 0.0
    return float(stats.ttest_rel(active[keep], reference[keep]).pvalue)


def independent_t_pvalue(active: np.ndarray, reference: np.ndarray) -> float:
    active = np.asarray(active, dtype=float)
    reference = np.asarray(reference, dtype=float)
    active = active[np.isfinite(active)]
    reference = reference[np.isfinite(reference)]
    if active.size < 2 or reference.size < 2:
        return np.nan
    return float(stats.ttest_ind(active, reference, equal_var=False).pvalue)


def plot_off_subject_diagnostic(
    trials: pd.DataFrame,
    out_base: Path,
    *,
    alpha: float,
) -> None:
    off = trials.loc[trials["session"].eq(1) & trials["medication"].eq("OFF") & np.isfinite(trials["inverse_rt_z"])].copy()
    if off.empty:
        return

    condition_codes = condition_code_order(off)
    label_lookup = (
        off.drop_duplicates("condition_code").set_index("condition_code")["condition_label"].astype(str).to_dict()
    )
    labels = [condition_display_label(label_lookup.get(code, code), code) for code in condition_codes]
    x = np.arange(len(condition_codes))

    subject_condition = (
        off.groupby(["subject", "condition_code"], observed=True)
        .agg(mean_inverse_rt_z=("inverse_rt_z", "mean"), sem_inverse_rt_z=("inverse_rt_z", "sem"))
        .reset_index()
    )
    subject_wide = subject_condition.pivot(index="subject", columns="condition_code", values="mean_inverse_rt_z")

    raw_p_by_code: dict[str, float] = {}
    q_by_code: dict[str, float] = {}
    if "gvs-01" in subject_wide.columns:
        raw_ps = []
        active_codes = [code for code in condition_codes if code != "gvs-01"]
        for code in active_codes:
            raw_ps.append(paired_t_pvalue(subject_wide[code].to_numpy(), subject_wide["gvs-01"].to_numpy()))
        finite = np.isfinite(raw_ps)
        q_values = np.full(len(raw_ps), np.nan, dtype=float)
        if np.any(finite):
            _, q_values[finite] = fdrcorrection(np.asarray(raw_ps, dtype=float)[finite], alpha=alpha)
        raw_p_by_code.update(dict(zip(active_codes, raw_ps, strict=True)))
        q_by_code.update(dict(zip(active_codes, q_values, strict=True)))

    subjects = sorted(off["subject"].unique(), key=lambda value: int(re.search(r"(\d+)$", value).group(1)))
    n_panels = len(subjects) + 2
    n_cols = 5
    n_rows = int(np.ceil(n_panels / n_cols))

    with plt.rc_context(STYLE):
        fig, axes = plt.subplots(n_rows, n_cols, figsize=(16.0, 3.0 * n_rows), constrained_layout=True)
        axes = np.asarray(axes).ravel()
        bar_colors = {"sham": "#FF0000", "active": "#0072BD", "significant": "#00FF00"}

        def draw_panel(
            ax: plt.Axes,
            means: pd.Series,
            sems: pd.Series,
            title: str,
            significant_codes: set[str],
            p_lookup: dict[str, float] | None = None,
        ) -> None:
            heights = np.asarray([means.get(code, np.nan) for code in condition_codes], dtype=float)
            errors = np.asarray([sems.get(code, np.nan) for code in condition_codes], dtype=float)
            colors = [
                bar_colors["sham"]
                if code == "gvs-01"
                else bar_colors["significant"]
                if code in significant_codes
                else bar_colors["active"]
                for code in condition_codes
            ]
            ax.bar(x, heights, yerr=errors, color=colors, edgecolor="black", linewidth=0.45, capsize=2)
            ax.axhline(0.0, color="black", linewidth=0.6)
            ax.set_title(title, fontsize=8, fontweight="normal", pad=2)
            ax.set_xticks(x)
            ax.set_xticklabels(labels, rotation=35, ha="right", fontsize=7)
            ax.tick_params(axis="y", labelsize=7)
            ax.grid(False)
            if p_lookup:
                y_min, y_max = ax.get_ylim()
                span = max(y_max - y_min, 1.0)
                for i, code in enumerate(condition_codes):
                    p_value = p_lookup.get(code, np.nan)
                    if np.isfinite(p_value) and p_value < alpha:
                        ax.text(i, y_max - 0.12 * span, f"{p_value:.3f}", ha="center", va="top", fontsize=7)

        for ax, subject in zip(axes, subjects, strict=False):
            subject_data = off.loc[off["subject"].eq(subject)]
            means = subject_data.groupby("condition_code", observed=True)["inverse_rt_z"].mean()
            sems = subject_data.groupby("condition_code", observed=True)["inverse_rt_z"].sem()
            sham_values = subject_data.loc[subject_data["condition_code"].eq("gvs-01"), "inverse_rt_z"].to_numpy()
            subject_ps = {
                code: independent_t_pvalue(
                    subject_data.loc[subject_data["condition_code"].eq(code), "inverse_rt_z"].to_numpy(),
                    sham_values,
                )
                for code in condition_codes
                if code != "gvs-01"
            }
            significant_codes = {code for code, p_value in subject_ps.items() if np.isfinite(p_value) and p_value < alpha}
            draw_panel(ax, means, sems, subject.replace("sub-pd", "PSPD"), significant_codes, subject_ps)

        all_subject_means = subject_condition.groupby("condition_code", observed=True)["mean_inverse_rt_z"].mean()
        all_subject_sems = subject_condition.groupby("condition_code", observed=True)["mean_inverse_rt_z"].sem()
        raw_significant = {code for code, p_value in raw_p_by_code.items() if np.isfinite(p_value) and p_value < alpha}
        fdr_significant = {code for code, q_value in q_by_code.items() if np.isfinite(q_value) and q_value < alpha}
        draw_panel(
            axes[len(subjects)],
            all_subject_means,
            all_subject_sems,
            "All Subjects/without FDR correction",
            raw_significant,
            raw_p_by_code,
        )
        draw_panel(
            axes[len(subjects) + 1],
            all_subject_means,
            all_subject_sems,
            "All Subjects/with FDR correction",
            fdr_significant,
            q_by_code,
        )

        for ax in axes[n_panels:]:
            ax.axis("off")
        fig.suptitle("1/RT, OFF", fontsize=12, fontweight="normal")
        fig.savefig(out_base.with_name(f"{out_base.name}_session1_off_subject_panels.png"), dpi=300, bbox_inches="tight")
        fig.savefig(out_base.with_name(f"{out_base.name}_session1_off_subject_panels.pdf"), bbox_inches="tight")
        plt.close(fig)


def plot_boxplots(
    trials: pd.DataFrame,
    anova: pd.DataFrame,
    pairwise: pd.DataFrame,
    out_base: Path,
    *,
    alpha: float,
) -> None:
    plot_data = trials.loc[np.isfinite(trials["rt_ms"])].copy()
    plot_data["normalized_rt"] = plot_data.groupby(["subject", "session"], group_keys=False)["rt_ms"].transform(
        zscore_series
    )
    plot_data = plot_data.loc[np.isfinite(plot_data["normalized_rt"])].copy()
    panels = (
        plot_data[["session", "medication"]]
        .drop_duplicates()
        .sort_values(["session", "medication"])
        .itertuples(index=False, name=None)
    )
    panels = list(panels)
    if len(panels) == 0:
        raise RuntimeError("No finite normalized RT values are available for plotting.")

    with plt.rc_context(STYLE):
        fig, axes = plt.subplots(1, len(panels), figsize=(6.15 * len(panels), 4.45), sharey=True)
        if len(panels) == 1:
            axes = [axes]

        positions = np.arange(len(CONDITION_ORDER))
        base_colors = {
            "sham": "#D1D5DB",
            "active": "#8AB6D6",
            "significant": "#D97706",
        }
        tick_fontsize = 13
        label_fontsize = 15
        panel_fontsize = 16
        box_linewidth = 1.8
        y_extents: list[tuple[float, float]] = []

        for panel_index, (ax, (session, medication)) in enumerate(zip(axes, panels, strict=True)):
            panel_label = chr(ord("A") + panel_index)
            subset = plot_data.loc[
                plot_data["session"].eq(session) & plot_data["medication"].eq(medication)
            ]
            values = [
                subset.loc[subset["condition_label"].eq(condition), "normalized_rt"].dropna().to_numpy(dtype=float)
                for condition in CONDITION_ORDER
            ]
            sig = pairwise.loc[
                pairwise["session"].eq(session)
                & pairwise["medication"].eq(medication)
                & pairwise["significant_raw_p"].astype(bool)
            ]
            significant_conditions = set(sig["condition_label"].astype(str))

            means = np.asarray([np.nanmean(item) if item.size else np.nan for item in values], dtype=float)
            sems = np.asarray(
                [stats.sem(item, nan_policy="omit") if np.count_nonzero(np.isfinite(item)) > 1 else np.nan for item in values],
                dtype=float,
            )
            colors = []
            for condition in CONDITION_ORDER:
                if condition == "sham":
                    colors.append(base_colors["sham"])
                elif condition in significant_conditions:
                    colors.append(base_colors["significant"])
                else:
                    colors.append(base_colors["active"])

            ax.bar(
                positions,
                means,
                yerr=sems,
                width=0.68,
                color=colors,
                edgecolor=colors,
                linewidth=box_linewidth,
                capsize=3,
                error_kw={"elinewidth": 1.0, "ecolor": "#374151", "capthick": 1.0},
            )

            ax.text(-0.075, 1.02, panel_label, transform=ax.transAxes, fontsize=panel_fontsize, fontweight="bold")
            ax.set_xticks(positions)
            code_lookup = (
                subset.drop_duplicates("condition_label")
                .set_index("condition_label")["condition_code"]
                .astype(str)
                .to_dict()
            )
            ax.set_xticklabels(
                [condition_display_label(condition, code_lookup.get(condition)) for condition in CONDITION_ORDER],
                rotation=35,
                ha="right",
                fontsize=tick_fontsize,
                fontweight="bold",
            )
            ax.tick_params(axis="y", labelsize=tick_fontsize)
            for tick_label in ax.get_yticklabels():
                tick_label.set_fontweight("bold")
            ax.axhline(0.0, color="#4B5563", linestyle="--", linewidth=0.9, zorder=0)
            lower = np.nanmin(means - np.nan_to_num(sems, nan=0.0))
            upper = np.nanmax(means + np.nan_to_num(sems, nan=0.0))
            y_extents.append((lower, upper))
            ax.set_xlabel("")
            ax.grid(False)

        lower = min(item[0] for item in y_extents)
        upper = max(item[1] for item in y_extents)
        padding = max(0.005, 0.03 * (upper - lower))
        axes[0].set_ylim(min(0.0, lower) - padding, max(0.0, upper) + padding)
        axes[0].set_ylabel("normalized RT", fontsize=label_fontsize, fontweight="bold")
        fig.subplots_adjust(wspace=0.10, left=0.07, right=0.995, bottom=0.24, top=0.95)
        fig.savefig(out_base.with_suffix(".png"), dpi=300, bbox_inches="tight")
        fig.savefig(out_base.with_suffix(".pdf"), bbox_inches="tight")
        plt.close(fig)


def write_summary(out_dir: Path, trials: pd.DataFrame, anova: pd.DataFrame, pairwise: pd.DataFrame) -> None:
    summary = {
        "n_trial_rows_loaded": int(len(trials)),
        "n_trial_rows_modelled": int(np.isfinite(trials["inverse_rt_z"]).sum()),
        "n_subjects": int(trials["subject"].nunique()),
        "conditions": CONDITION_ORDER,
        "analysis_metric": "1/RT z-scored within subject/session",
        "anova": anova.to_dict(orient="records"),
        "significant_gvs_vs_sham_fdr": pairwise.loc[pairwise["significant_fdr"].astype(bool)].to_dict(orient="records"),
    }
    (out_dir / "gvs_rt_lme_summary.json").write_text(json.dumps(summary, indent=2, default=json_default))


def json_default(value: Any) -> Any:
    if isinstance(value, (np.integer,)):
        return int(value)
    if isinstance(value, (np.floating,)):
        return float(value)
    if pd.isna(value):
        return None
    return str(value)


def main() -> None:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    if args.data_source == "metrics_mat":
        trials = load_trial_table_from_metrics(args.metrics_dir, args.input_kind, args.rt_column)
        mapping = "behav_metrics cell order: sham, GVS1, GVS2, ..., GVS8"
        source_path = args.metrics_dir
    elif args.use_inventory_trial_ranges:
        trials = load_trial_table(args.behaviour_dir, args.inventory, args.input_kind, args.rt_column)
        mapping = "inventory trial_start/trial_stop ranges"
        source_path = args.behaviour_dir
    else:
        trials = load_trial_table_from_stim_order(args.behaviour_dir, args.inventory, args.input_kind, args.rt_column)
        mapping = "stim_order full local RT blocks"
        source_path = args.behaviour_dir
    anova, pairwise, terms = fit_session_models(
        trials,
        include_run=not args.no_run_fixed_effect,
        alpha=args.alpha,
    )
    condition_summary = summarize_conditions(trials)

    trials.to_csv(args.out_dir / "gvs_rt_lme_trial_table.csv", index=False)
    condition_summary.to_csv(args.out_dir / "gvs_rt_lme_condition_summary.csv", index=False)
    anova.to_csv(args.out_dir / "gvs_rt_lme_anova_condition_effect.csv", index=False)
    pairwise.to_csv(args.out_dir / "gvs_rt_lme_pairwise_gvs_vs_sham.csv", index=False)
    terms.to_csv(args.out_dir / "gvs_rt_lme_model_terms.csv", index=False)
    plot_boxplots(trials, anova, pairwise, args.out_dir / "gvs_rt_lme_boxplot", alpha=args.alpha)
    plot_off_subject_diagnostic(trials, args.out_dir / "gvs_rt_lme_boxplot", alpha=args.alpha)
    write_summary(args.out_dir, trials, anova, pairwise)

    print(f"Loaded {len(trials)} trial rows from {source_path}")
    print(f"Trial mapping: {mapping}")
    print(f"Modelled {int(np.isfinite(trials['inverse_rt_z']).sum())} finite z-scored 1/RT trials")
    for _, row in anova.iterrows():
        sig_raw = pairwise.loc[
            pairwise["session"].eq(row["session"])
            & pairwise["medication"].eq(row["medication"])
            & pairwise["significant_raw_p"].astype(bool)
        ]["condition_label"].tolist()
        sig_fdr = pairwise.loc[
            pairwise["session"].eq(row["session"])
            & pairwise["medication"].eq(row["medication"])
            & pairwise["significant_fdr"].astype(bool)
        ]["condition_label"].tolist()
        sig_raw_text = ", ".join(sig_raw) if sig_raw else "none"
        sig_fdr_text = ", ".join(sig_fdr) if sig_fdr else "none"
        print(
            f"Session {int(row['session'])} {row['medication']}: "
            f"ANOVA chi2({int(row['df'])}) = {row['chi2']:.3f}, "
            f"{format_p(float(row['p_value']))}; "
            f"raw p<{args.alpha:g}: {sig_raw_text}; FDR q<{args.alpha:g}: {sig_fdr_text}"
        )
    print(f"Wrote outputs to {args.out_dir}")


if __name__ == "__main__":
    main()
