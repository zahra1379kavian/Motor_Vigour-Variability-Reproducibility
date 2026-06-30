#!/usr/bin/env python3
"""Refresh curated Figure 8 panels from packaged result tables and exports."""

from __future__ import annotations

import argparse
import shutil
import warnings
from pathlib import Path

warnings.filterwarnings("ignore", message="Unable to import Axes3D.*")

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_RESULTS_DIR = ROOT / "results" / "main" / "figure_08_gvs_effects"
DEFAULT_FIGURES_DIR = ROOT / "figures" / "main"

CONDITION_ORDER = ["sham"] + [f"GVS{i}" for i in range(1, 9)]
SHAM_CODE = "gvs-01"
FEATURE = "early_late_change"
FEATURE_LABEL = "Late-Early Response Difference (GVS-Sham)"
SIGNIFICANCE_Q_THRESHOLD = 0.06

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
    parser.add_argument("--results-dir", type=Path, default=DEFAULT_RESULTS_DIR)
    parser.add_argument("--figures-dir", type=Path, default=DEFAULT_FIGURES_DIR)
    parser.add_argument(
        "--replot-feature-panels",
        action="store_true",
        help=(
            "Regenerate Figure 8B/8C from packaged CSV tables. By default the "
            "exact curated 8B/8C source exports are copied into figures/main."
        ),
    )
    return parser.parse_args()


def zscore_series(values: pd.Series) -> pd.Series:
    numeric = pd.to_numeric(values, errors="coerce")
    sd = numeric.std(ddof=0)
    if not np.isfinite(sd) or sd == 0:
        return pd.Series(np.nan, index=values.index)
    return (numeric - numeric.mean()) / sd


def condition_display_label(condition_label: str, condition_code: str | None = None) -> str:
    label = str(condition_label)
    if label.lower() == "sham":
        return "Sham"
    if label.upper().startswith("GVS"):
        return label.upper()
    if condition_code:
        code = str(condition_code).lower().replace("_", "-")
        if code == SHAM_CODE:
            return "Sham"
        if code.startswith("gvs-"):
            try:
                return f"GVS{int(code.split('-', 1)[1]) - 1}"
            except ValueError:
                return label
    return label


def gvs_display_label(gvs_code: str) -> str:
    code = str(gvs_code).strip().lower().replace("_", "-").replace(" ", "-")
    if code in {SHAM_CODE, "sham"}:
        return "sham"
    if not code.startswith("gvs-"):
        return str(gvs_code)
    try:
        index = int(code.split("-", 1)[1])
    except ValueError:
        return str(gvs_code)
    return f"gvs{index - 1}" if index > 1 else "sham"


def plot_rt_panel(results_dir: Path, figures_dir: Path) -> None:
    rt_dir = results_dir / "rt_lme"
    trials = pd.read_csv(rt_dir / "gvs_rt_lme_trial_table.csv")
    pairwise = pd.read_csv(rt_dir / "gvs_rt_lme_pairwise_gvs_vs_sham.csv")

    plot_data = trials.loc[np.isfinite(trials["rt_ms"])].copy()
    plot_data["normalized_rt"] = plot_data.groupby(["subject", "session"], group_keys=False)["rt_ms"].transform(
        zscore_series
    )
    plot_data = plot_data.loc[np.isfinite(plot_data["normalized_rt"])].copy()
    panels = list(
        plot_data[["session", "medication"]]
        .drop_duplicates()
        .sort_values(["session", "medication"])
        .itertuples(index=False, name=None)
    )
    if not panels:
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
        y_extents: list[tuple[float, float]] = []

        for panel_index, (ax, (session, medication)) in enumerate(zip(axes, panels, strict=True)):
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
                [
                    stats.sem(item, nan_policy="omit") if np.count_nonzero(np.isfinite(item)) > 1 else np.nan
                    for item in values
                ],
                dtype=float,
            )
            colors = [
                base_colors["sham"]
                if condition == "sham"
                else base_colors["significant"]
                if condition in significant_conditions
                else base_colors["active"]
                for condition in CONDITION_ORDER
            ]

            ax.bar(
                positions,
                means,
                yerr=sems,
                width=0.68,
                color=colors,
                edgecolor=colors,
                linewidth=1.8,
                capsize=3,
                error_kw={"elinewidth": 1.0, "ecolor": "#374151", "capthick": 1.0},
            )
            ax.text(-0.075, 1.02, chr(ord("A") + panel_index), transform=ax.transAxes, fontsize=16, fontweight="bold")
            code_lookup = (
                subset.drop_duplicates("condition_label")
                .set_index("condition_label")["condition_code"]
                .astype(str)
                .to_dict()
            )
            ax.set_xticks(positions)
            ax.set_xticklabels(
                [condition_display_label(condition, code_lookup.get(condition)) for condition in CONDITION_ORDER],
                rotation=35,
                ha="right",
                fontsize=13,
                fontweight="bold",
            )
            ax.tick_params(axis="y", labelsize=13)
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
        axes[0].set_ylabel("normalized RT", fontsize=15, fontweight="bold")
        fig.subplots_adjust(wspace=0.10, left=0.07, right=0.995, bottom=0.24, top=0.95)
        figures_dir.mkdir(parents=True, exist_ok=True)
        fig.savefig(figures_dir / "figure_08a_gvs_rt_lme_boxplot.png", dpi=300, bbox_inches="tight")
        fig.savefig(figures_dir / "figure_08a_gvs_rt_lme_boxplot.pdf", bbox_inches="tight")
        plt.close(fig)


def make_session_subject_pairs(run_features: pd.DataFrame, session: int, medication: str) -> pd.DataFrame:
    session_runs = run_features.loc[
        run_features["session"].astype(int).eq(int(session))
        & run_features["medication"].astype(str).str.upper().eq(medication.upper())
    ].copy()
    if session_runs.empty:
        raise ValueError(f"No run-level features found for session {session} / {medication}")

    index_cols = ["subject", "session", "medication", "run"]
    sham_runs = session_runs.loc[session_runs["gvs_code"].eq(SHAM_CODE)].set_index(index_cols)
    active_codes = sorted(code for code in session_runs["gvs_code"].unique() if code != SHAM_CODE)
    rows = []
    for active_code in active_codes:
        active_runs = session_runs.loc[session_runs["gvs_code"].eq(active_code)].set_index(index_cols)
        paired_runs = sham_runs.join(active_runs, how="inner", lsuffix="_sham", rsuffix="_active")
        if paired_runs.empty:
            continue
        per_run = pd.DataFrame(
            {
                "subject": paired_runs.index.get_level_values("subject"),
                "sham_value": paired_runs[f"{FEATURE}_sham"].to_numpy(dtype=np.float64),
                "active_value": paired_runs[f"{FEATURE}_active"].to_numpy(dtype=np.float64),
            }
        )
        per_run["delta_active_minus_sham"] = per_run["active_value"] - per_run["sham_value"]
        per_subject = (
            per_run.replace([np.inf, -np.inf], np.nan)
            .dropna(subset=["sham_value", "active_value", "delta_active_minus_sham"])
            .groupby("subject", as_index=False)
            .agg(
                sham_value=("sham_value", "mean"),
                active_value=("active_value", "mean"),
                delta_active_minus_sham=("delta_active_minus_sham", "mean"),
                n_runs=("delta_active_minus_sham", "size"),
            )
        )
        per_subject.insert(0, "feature", FEATURE)
        per_subject.insert(0, "active_gvs", active_code)
        rows.append(per_subject)

    if not rows:
        raise RuntimeError(f"No active-vs-sham subject pairs could be built against {SHAM_CODE}")
    return pd.concat(rows, ignore_index=True).sort_values(["active_gvs", "subject"]).reset_index(drop=True)


def padded_limits(values: list[np.ndarray]) -> tuple[float, float]:
    finite = np.concatenate([np.asarray(v, dtype=float) for v in values])
    finite = finite[np.isfinite(finite)]
    if finite.size == 0:
        return -1.0, 1.0
    lo = min(0.0, float(np.nanmin(finite)))
    hi = max(0.0, float(np.nanmax(finite)))
    span = hi - lo
    if span == 0.0:
        span = max(abs(lo), abs(hi), 1.0)
    return lo - 0.08 * span, hi + 0.08 * span


def plot_feature_panel(input_dir: Path, output_png: Path, output_pdf: Path) -> None:
    run_features = pd.read_csv(input_dir / "gvs_run_signal_features.csv")
    sessions = [
        (1, "OFF", "session1_off"),
        (2, "ON", "session2_on"),
    ]
    panel_pairs = []
    stats_lookup = {}
    for session, medication, suffix in sessions:
        subject_pairs = make_session_subject_pairs(run_features, session, medication)
        panel_pairs.append(subject_pairs)
        stats_path = input_dir / f"gvs_vs_sham_signal_feature_stats_{suffix}.csv"
        stats_df = pd.read_csv(stats_path) if stats_path.exists() else pd.DataFrame()
        if not stats_df.empty:
            for row in stats_df.itertuples(index=False):
                if getattr(row, "feature") == FEATURE:
                    stats_lookup[(session, getattr(row, "active_gvs"))] = float(getattr(row, "q_perm_fdr"))

    active_codes = sorted(panel_pairs[0]["active_gvs"].unique())
    x = np.arange(len(active_codes))
    rng = np.random.default_rng(0)

    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Liberation Sans", "Arial", "DejaVu Sans"],
            "font.size": 12,
            "font.weight": "bold",
            "axes.labelweight": "bold",
            "axes.titleweight": "bold",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )
    fig, axes = plt.subplots(1, 2, figsize=(13.2, 3.7), constrained_layout=True)
    fig.set_constrained_layout_pads(w_pad=0.02, h_pad=0.02, wspace=0.04, hspace=0.04)

    for ax, (subject_pairs, (session, _medication, _suffix)) in zip(axes, zip(panel_pairs, sessions), strict=True):
        groups = [
            subject_pairs.loc[
                subject_pairs["active_gvs"].eq(code) & subject_pairs["feature"].eq(FEATURE),
                "delta_active_minus_sham",
            ]
            .dropna()
            .to_numpy(dtype=np.float64)
            for code in active_codes
        ]
        significant_groups = [
            np.isfinite(stats_lookup.get((session, code), np.nan))
            and float(stats_lookup.get((session, code), np.nan)) < SIGNIFICANCE_Q_THRESHOLD
            for code in active_codes
        ]
        violins = ax.violinplot(
            groups,
            positions=x,
            widths=0.74,
            showmeans=False,
            showmedians=False,
            showextrema=False,
        )
        for body, is_significant in zip(violins["bodies"], significant_groups):
            body.set_facecolor("#d98c4a" if is_significant else "#8fb3c8")
            body.set_edgecolor("#8f4d1f" if is_significant else "#466f86")
            body.set_alpha(0.72)
            body.set_linewidth(1.0 if is_significant else 0.8)

        box = ax.boxplot(
            groups,
            positions=x,
            widths=0.28,
            patch_artist=True,
            showfliers=False,
            medianprops={"color": "#111111", "linewidth": 1.2},
            whiskerprops={"color": "#333333", "linewidth": 0.8},
            capprops={"color": "#333333", "linewidth": 0.8},
        )
        for patch, is_significant in zip(box["boxes"], significant_groups):
            patch.set_facecolor("#f4c08d" if is_significant else "#f6f6f2")
            patch.set_edgecolor("#8f4d1f" if is_significant else "#333333")
            patch.set_linewidth(1.0 if is_significant else 0.8)

        for xpos, values in zip(x, groups):
            jitter = rng.uniform(-0.12, 0.12, size=values.size)
            ax.scatter(
                np.full(values.size, xpos) + jitter,
                values,
                s=19,
                color="#222222",
                alpha=0.68,
                linewidths=0,
                zorder=3,
            )

        y_min, y_max = padded_limits(groups)
        ax.set_ylim(y_min, y_max)
        ax.axhline(0.0, color="#222222", linewidth=0.8)
        ax.set_xticks(x)
        ax.set_xlim(x[0] - 0.5, x[-1] + 0.5)
        ax.set_xticklabels(
            [gvs_display_label(code) for code in active_codes],
            rotation=0,
            fontsize=12,
            fontweight="bold",
        )
        ax.tick_params(axis="both", labelsize=12)
        for label in ax.get_yticklabels():
            label.set_fontweight("bold")
        ax.set_ylabel(FEATURE_LABEL, fontsize=13, fontweight="bold")
        ax.grid(False)
        ax.spines[["top", "right"]].set_visible(False)

    output_png.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_png, dpi=300)
    fig.savefig(output_pdf)
    plt.close(fig)


def copy_curated_feature_panel(input_dir: Path, output_png: Path, output_pdf: Path) -> None:
    output_png.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(input_dir / "gvs_vs_sham_feature_delta_first_subplot_off_on.png", output_png)
    shutil.copy2(input_dir / "gvs_vs_sham_feature_delta_first_subplot_off_on.pdf", output_pdf)


def main() -> int:
    args = parse_args()
    plot_rt_panel(args.results_dir, args.figures_dir)
    feature_panel_writer = plot_feature_panel if args.replot_feature_panels else copy_curated_feature_panel
    feature_panel_writer(
        args.results_dir / "vigour_network_feature_delta",
        args.figures_dir / "figure_08b_gvs_vigour_network_feature_delta.png",
        args.figures_dir / "figure_08b_gvs_vigour_network_feature_delta.pdf",
    )
    feature_panel_writer(
        args.results_dir / "task_activation_feature_delta",
        args.figures_dir / "figure_08c_gvs_task_activation_feature_delta.png",
        args.figures_dir / "figure_08c_gvs_task_activation_feature_delta.pdf",
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
