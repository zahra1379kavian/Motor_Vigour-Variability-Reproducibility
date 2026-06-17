#!/usr/bin/env python3
"""Compare intra-vs-between FC medication effects with run-to-run baselines.

This is a companion to ``med_effects.py`` and the existing
``intra_vs_between_fc_medication_change.png`` figure.  The original result uses
concatenated runs within each session and compares ON minus OFF.  Here we
estimate the same intra-ROI voxel-pair FC and between-ROI FC separately for run
1 and run 2, then use run2 - run1 inside each medication state as a within-state
baseline for comparison with the current ON - OFF result.
"""


import argparse
import json
import re
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import nibabel as nib
import numpy as np
import pandas as pd
from matplotlib.text import Text

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
for path in (ROOT, HERE):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

import med_effects as M  # noqa: E402

DEFAULT_OUT = ROOT / "results" / "supplementary" / "figure_12_run_baseline_fc" / "vigour_network"
RUN_RE = re.compile(r"_run-(?P<run>\d+)")
RUN_1 = "1"
RUN_2 = "2"


def build_parser():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--weight-map", type=Path, default=M.DEFAULT_WEIGHT_MAP)
    parser.add_argument("--roi-definition-figure", type=Path, default=M.DEFAULT_ROI_FIGURE)
    parser.add_argument("--roi-region-table", type=Path, default=None)
    parser.add_argument("--roi-percentile", type=float, default=M.REFERENCE_THRESHOLD)
    parser.add_argument("--min-report-voxels", type=int, default=M.DEFAULT_MIN_REPORT_VOXELS)
    parser.add_argument("--session-manifest", type=Path, default=M.DEFAULT_SESSION_MANIFEST)
    parser.add_argument("--beta-root", type=Path, default=M.DEFAULT_BETA_ROOT)
    parser.add_argument("--subjects", nargs="+", default=None)
    parser.add_argument("--complete-subjects-only", action="store_true")
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--voxel-selection", choices=M.VOXEL_SELECTION_MODES, default=M.VOXEL_SELECTION_WEIGHTED_VIGOUR)
    hemisphere_group = parser.add_mutually_exclusive_group()
    hemisphere_group.add_argument("--split-hemispheres", dest="split_hemispheres", action="store_true", default=True)
    hemisphere_group.add_argument("--no-split-hemispheres", dest="split_hemispheres", action="store_false")
    parser.add_argument("--exclude-rois", nargs="*", default=())
    parser.add_argument("--min-lateralized-voxels", type=int, default=1)
    parser.add_argument("--aal-version", default=M.DEFAULT_AAL_VERSION)
    parser.add_argument("--atlas-cache-dir", type=Path, default=M.DEFAULT_ATLAS_CACHE_DIR)
    parser.add_argument("--intra-between-fc-metric", choices=M.INTRA_BETWEEN_FC_METRICS, default=M.INTRA_BETWEEN_FC_METRIC)
    parser.add_argument("--mi-quantile-bins", type=int, default=M.DEFAULT_MI_QUANTILE_BINS)
    parser.add_argument("--random-state", type=int, default=0)
    return parser


def _run_from_beta_path(path):
    match = M.BETA_FILE_RE.match(path.name)
    if match is not None:
        return str(int(match.group("run")))
    match = RUN_RE.search(path.name)
    if match is not None:
        return str(int(match.group("run")))
    raise ValueError(f"Could not parse run number from {path}")


def _single_run_specs(specs):
    run_specs = []
    missing = []
    for spec in specs:
        if not spec.beta_paths:
            missing.append(spec.label)
            continue
        for beta_path in spec.beta_paths:
            run = _run_from_beta_path(beta_path)
            label = f"{spec.label}_run-{run}"
            run_spec = M.SessionSpec(
                label=label,
                subject=spec.subject,
                session=spec.session,
                state=spec.state,
                bold_path=None,
                timeseries_path=None,
                beta_paths=(beta_path,),
            )
            run_specs.append((run_spec, run))
    if missing:
        raise ValueError(
            "Run-baseline analysis requires beta_path inputs; missing beta paths for: "
            + ", ".join(missing)
        )
    return run_specs


def _prepare_rois(args):
    if args.roi_region_table is None:
        args.roi_region_table = M._default_region_table_for(args.roi_definition_figure)
    missing = M._missing_inputs(args)
    if missing:
        raise ValueError("Missing required inputs:\n- " + "\n- ".join(missing))

    weight_img = nib.load(str(args.weight_map))
    weight_values = np.asarray(weight_img.get_fdata(), dtype=np.float64)
    groups, _, roi_names, min_roi_voxels = M._analysis_roi_setup(args, weight_img)
    rois, roi_threshold, weighted_rois = M._build_analysis_rois(
        weight_values=weight_values,
        roi_names=roi_names,
        groups=groups,
        roi_percentile=args.roi_percentile,
        min_report_voxels=args.min_report_voxels,
        min_roi_voxels=min_roi_voxels,
        voxel_selection=args.voxel_selection,
        random_state=args.random_state,
    )
    return weight_img, rois, weighted_rois, roi_threshold


def _compute_run_level_values(args, weight_img, rois):
    specs = M._load_session_specs(args)
    run_specs = _single_run_specs(specs)
    session_rows = []
    roi_rows = []
    for run_spec, run in run_specs:
        roi_ts = M._load_session_timeseries(run_spec, weight_img, rois)
        cleaned = M._clean_timeseries(roi_ts)
        row = {
            "label": run_spec.label,
            "subject": run_spec.subject,
            "session": run_spec.session,
            "state": run_spec.state,
            "run": run,
            "connectivity_metric": M.INTRA_BETWEEN_FC_METRIC,
        }
        row.update(M._between_roi_fc_summary(cleaned, mi_quantile_bins=args.mi_quantile_bins))
        voxel_ts = M._load_session_voxel_timeseries(run_spec, weight_img, rois)
        run_roi_rows, intra_summary = M._intra_roi_fc_values(
            run_spec,
            voxel_ts,
            rois,
            mi_quantile_bins=args.mi_quantile_bins,
        )
        for roi_row in run_roi_rows:
            roi_row["run"] = run
        row.update(intra_summary)
        session_rows.append(row)
        roi_rows.extend(run_roi_rows)
        print(f"Computed {run_spec.label}", flush=True)
    values = pd.DataFrame(session_rows).sort_values(["subject", "session", "run"]).reset_index(drop=True)
    roi_values = pd.DataFrame(roi_rows).sort_values(["subject", "session", "run", "roi"]).reset_index(drop=True)
    return values, roi_values


def _metric_columns(values):
    if "within_roi_mean_mi" in values.columns or "within_roi_delta_mi_on_minus_off" in values.columns:
        return {
            "within": "within_roi_mean_mi",
            "between": "between_roi_mean_mi",
            "current_within_delta": "within_roi_delta_mi_on_minus_off",
            "current_between_delta": "between_roi_delta_mi_on_minus_off",
            "current_contrast_delta": "within_minus_between_delta_mi",
            "run_within_delta": "within_roi_delta_mi_run2_minus_run1",
            "run_between_delta": "between_roi_delta_mi_run2_minus_run1",
            "run_contrast_delta": "within_minus_between_delta_mi_run2_minus_run1",
            "plot_ylabel": "FC change",
            "metric_suffix": "mi",
        }
    return {
        "within": "within_roi_mean_z",
        "between": "between_roi_mean_z",
        "current_within_delta": "within_roi_delta_z_on_minus_off",
        "current_between_delta": "between_roi_delta_z_on_minus_off",
        "current_contrast_delta": "within_minus_between_delta_z",
        "run_within_delta": "within_roi_delta_z_run2_minus_run1",
        "run_between_delta": "between_roi_delta_z_run2_minus_run1",
        "run_contrast_delta": "within_minus_between_delta_z_run2_minus_run1",
        "plot_ylabel": "FC change",
        "metric_suffix": "z",
    }


def _complete_run_deltas(run_values):
    cols = _metric_columns(run_values)
    rows = []
    for (subject, state), group in run_values.groupby(["subject", "state"], sort=True):
        run1 = group.loc[group["run"].astype(str) == RUN_1]
        run2 = group.loc[group["run"].astype(str) == RUN_2]
        if run1.shape[0] != 1 or run2.shape[0] != 1:
            continue
        run1 = run1.iloc[0]
        run2 = run2.iloc[0]
        within_delta = float(run2[cols["within"]] - run1[cols["within"]])
        between_delta = float(run2[cols["between"]] - run1[cols["between"]])
        row = {
            "subject": subject,
            "state": state,
            "session": run1["session"],
            "run1_label": run1["label"],
            "run2_label": run2["label"],
            cols["run_within_delta"]: within_delta,
            cols["run_between_delta"]: between_delta,
            cols["run_contrast_delta"]: within_delta - between_delta,
        }
        if cols["metric_suffix"] == "z":
            row.update(
                {
                    "within_roi_run1_z": float(run1[cols["within"]]),
                    "within_roi_run2_z": float(run2[cols["within"]]),
                    "within_roi_run1_r": float(np.tanh(run1[cols["within"]])),
                    "within_roi_run2_r": float(np.tanh(run2[cols["within"]])),
                    "within_roi_delta_r_run2_minus_run1": float(np.tanh(run2[cols["within"]]) - np.tanh(run1[cols["within"]])),
                    "between_roi_run1_z": float(run1[cols["between"]]),
                    "between_roi_run2_z": float(run2[cols["between"]]),
                    "between_roi_run1_r": float(np.tanh(run1[cols["between"]])),
                    "between_roi_run2_r": float(np.tanh(run2[cols["between"]])),
                    "between_roi_delta_r_run2_minus_run1": float(np.tanh(run2[cols["between"]]) - np.tanh(run1[cols["between"]])),
                }
            )
        else:
            row.update(
                {
                    "within_roi_run1_mi": float(run1[cols["within"]]),
                    "within_roi_run2_mi": float(run2[cols["within"]]),
                    "between_roi_run1_mi": float(run1[cols["between"]]),
                    "between_roi_run2_mi": float(run2[cols["between"]]),
                }
            )
        rows.append(row)
    return pd.DataFrame(rows).sort_values(["subject", "state"]).reset_index(drop=True)


def _run_level_medication_effect(group, cols):
    values = {}
    for state in ("off", "on"):
        state_values = group.loc[group["state"] == state]
        for run in (RUN_1, RUN_2):
            row = state_values.loc[state_values["run"].astype(str) == run]
            if row.shape[0] != 1:
                return None
            row = row.iloc[0]
            within = float(row[cols["within"]])
            between = float(row[cols["between"]])
            values[(state, run)] = (within, between, within - between)
    run1 = np.asarray(values[("on", RUN_1)], dtype=np.float64) - np.asarray(values[("off", RUN_1)], dtype=np.float64)
    run2 = np.asarray(values[("on", RUN_2)], dtype=np.float64) - np.asarray(values[("off", RUN_2)], dtype=np.float64)
    effect = 0.5 * (np.abs(run1) + np.abs(run2))
    return (float(effect[0]), float(effect[1]), float(effect[2]))


def _comparison_subject_values(run_values, run_deltas):
    cols = _metric_columns(run_values)
    rows = []
    for subject, group in run_deltas.groupby("subject", sort=True):
        off = group.loc[group["state"] == "off"]
        on = group.loc[group["state"] == "on"]
        if off.shape[0] != 1 or on.shape[0] != 1:
            continue
        off = off.iloc[0]
        on = on.iloc[0]
        medication_effect = _run_level_medication_effect(run_values.loc[run_values["subject"] == subject], cols)
        if medication_effect is None:
            continue
        within_baseline = 0.5 * (abs(float(off[cols["run_within_delta"]])) + abs(float(on[cols["run_within_delta"]])))
        between_baseline = 0.5 * (abs(float(off[cols["run_between_delta"]])) + abs(float(on[cols["run_between_delta"]])))
        contrast_baseline = 0.5 * (abs(float(off[cols["run_contrast_delta"]])) + abs(float(on[cols["run_contrast_delta"]])))
        source_specs = [
            (
                "run_to_run_baseline_abs",
                "Run-to-run baseline",
                within_baseline,
                between_baseline,
                contrast_baseline,
            ),
            (
                "run_level_medication_effect",
                "Run-level medication effect",
                medication_effect[0],
                medication_effect[1],
                medication_effect[2],
            ),
        ]
        for source_key, source_label, within_value, between_value, contrast_value in source_specs:
            rows.append(
                {
                    "subject": subject,
                    "source": source_key,
                    "source_label": source_label,
                    "within_roi_change": float(within_value),
                    "between_roi_change": float(between_value),
                    "within_minus_between_change": float(contrast_value),
                }
            )
    return pd.DataFrame(rows)


def _source_stats(comparison_values):
    components = [
        ("within_roi_change", "Intra-ROI FC change"),
        ("between_roi_change", "Between-ROI FC change"),
        ("within_minus_between_change", "Intra-minus-between FC change"),
    ]
    rows = []
    for source, group in comparison_values.groupby("source", sort=False):
        source_label = group["source_label"].iloc[0]
        for component, description in components:
            stats = M._single_subject_level_test(group, component)
            row = {
                "analysis_type": "source_vs_zero",
                "source": source,
                "source_label": source_label,
                "baseline_source": "",
                "component": component,
                "description": description,
            }
            row.update(stats)
            rows.append(row)
    return pd.DataFrame(rows)


def _paired_comparison_stats(comparison_values):
    components = [
        ("within_roi_change", "Medication-effect magnitude minus baseline for intra-ROI FC"),
        ("between_roi_change", "Medication-effect magnitude minus baseline for between-ROI FC"),
        ("within_minus_between_change", "Medication-effect magnitude minus baseline for intra-minus-between FC"),
    ]
    rows = []
    for component, description in components:
        wide = comparison_values.pivot(index="subject", columns="source", values=component)
        needed = wide.loc[:, ["run_to_run_baseline_abs", "run_level_medication_effect"]].dropna()
        values = pd.DataFrame(
            {
                "subject": needed.index,
                "medication_minus_baseline": needed["run_level_medication_effect"].to_numpy(dtype=np.float64)
                - needed["run_to_run_baseline_abs"].to_numpy(dtype=np.float64),
            }
        )
        stats = M._single_subject_level_test(values, "medication_minus_baseline")
        row = {
            "analysis_type": "medication_minus_baseline",
            "source": "run_level_medication_effect",
            "source_label": "Run-level medication effect",
            "baseline_source": "run_to_run_baseline_abs",
            "component": component,
            "description": description,
        }
        row.update(stats)
        rows.append(row)
    return pd.DataFrame(rows)


def _format_p(value):
    if not np.isfinite(float(value)):
        return "n/a"
    if float(value) < 0.001:
        return "<0.001"
    return f"{float(value):.3f}"


def _plot_comparison(comparison_values, summary, out_dir, ylabel):
    components = [
        ("within_roi_change", "Intra-ROI"),
        ("between_roi_change", "Between-ROI"),
        ("within_minus_between_change", "Intra - Between"),
    ]
    sources = [
        ("run_to_run_baseline_abs", "Between\nRuns", "#555555"),
        ("run_level_medication_effect", "Between\nSessions", "#D62728"),
    ]
    x_positions = [0.0, 0.58]
    fig, axes = plt.subplots(1, 3, figsize=(9.2, 3.6), sharey=False)
    rng = np.random.default_rng(0)
    for ax, (component, title) in zip(axes, components):
        wide = comparison_values.pivot(index="subject", columns="source", values=component)
        for _, row in wide.iterrows():
            if np.isfinite(row.get("run_to_run_baseline_abs", np.nan)) and np.isfinite(row.get("run_level_medication_effect", np.nan)):
                ax.plot(
                    x_positions,
                    [row["run_to_run_baseline_abs"], row["run_level_medication_effect"]],
                    color="#b3b3b3",
                    linewidth=0.8,
                    alpha=0.5,
                    zorder=1,
                )
        finite_values = []
        summary_bounds = []
        for x_value, (source, _, color) in zip(x_positions, sources):
            vals = comparison_values.loc[comparison_values["source"] == source, component].to_numpy(dtype=np.float64)
            vals = vals[np.isfinite(vals)]
            finite_values.extend(vals.tolist())
            ax.scatter(
                rng.normal(x_value, 0.018, vals.size),
                vals,
                s=42,
                color=color,
                edgecolor="white",
                linewidth=0.45,
                alpha=0.9,
                zorder=2,
            )
            stat = summary.loc[
                (summary["analysis_type"] == "source_vs_zero")
                & (summary["source"] == source)
                & (summary["component"] == component)
            ]
            if not stat.empty:
                stat = stat.iloc[0]
                mean_value = float(stat["mean"])
                ci_low = float(stat.get("ci95_low", np.nan))
                ci_high = float(stat.get("ci95_high", np.nan))
                if np.isfinite(ci_low) and np.isfinite(ci_high):
                    ax.errorbar(
                        x_value,
                        mean_value,
                        yerr=[[mean_value - ci_low], [ci_high - mean_value]],
                        fmt="o",
                        color="#111111",
                        ecolor="#111111",
                        elinewidth=2.0,
                        capsize=4.5,
                        markersize=6.2,
                        zorder=3,
                    )
                    summary_bounds.extend([ci_low, ci_high])
                else:
                    ax.scatter([x_value], [mean_value], s=54, color="#111111", zorder=3)
        comparison = summary.loc[
            (summary["analysis_type"] == "medication_minus_baseline")
            & (summary["baseline_source"] == "run_to_run_baseline_abs")
            & (summary["component"] == component)
        ]
        y_all = np.asarray(finite_values + summary_bounds, dtype=np.float64)
        if y_all.size:
            y_min = float(np.nanmin(y_all))
            y_max = float(np.nanmax(y_all))
        else:
            y_min, y_max = -1.0, 1.0
        y_range = y_max - y_min
        if not np.isfinite(y_range) or y_range <= 0:
            y_range = 1.0
        label_y = y_max + 0.12 * y_range
        ax.plot(x_positions, [label_y, label_y], color="#333333", linewidth=1.0, clip_on=False)
        p_value = float(comparison["paired_t_p_value_two_sided"].iloc[0]) if not comparison.empty else np.nan
        ax.text(float(np.mean(x_positions)), label_y + 0.025 * y_range, f"paired contrast p = {_format_p(p_value)}", ha="center", va="bottom", fontsize=11.0)
        ax.set_ylim(y_min - 0.12 * y_range, label_y + 0.22 * y_range)
        ax.axhline(0.0, color="#666666", linestyle="--", linewidth=0.9, alpha=0.7, zorder=0)
        ax.set_xticks(x_positions)
        ax.set_xticklabels([label for _, label, _ in sources])
        ax.set_xlim(-0.18, 0.76)
        ax.tick_params(labelsize=12)
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.set_axisbelow(True)
    axes[0].set_ylabel(ylabel)
    for label, ax in zip(("A", "B", "C"), axes):
        ax.text(-0.16, 1.04, label, transform=ax.transAxes, fontsize=20, fontweight="bold", ha="left", va="bottom")
    M._apply_paper_typography(fig, axes)
    for ax in axes:
        ax.tick_params(axis="both", labelsize=12)
        ax.yaxis.label.set_fontsize(15)
        ax.yaxis.label.set_fontweight("bold")
        for tick_label in ax.get_xticklabels():
            tick_label.set_fontsize(12)
            tick_label.set_fontweight("bold")
        for tick_label in ax.get_yticklabels():
            tick_label.set_fontsize(12)
            tick_label.set_fontweight("bold")
    for text in fig.findobj(match=Text):
        text.set_fontweight("bold")
    fig.tight_layout(w_pad=0.6)
    out_dir.mkdir(parents=True, exist_ok=True)
    png_path = out_dir / "intra_vs_between_fc_run_baseline_comparison.png"
    with plt.rc_context({"pdf.fonttype": 42, "ps.fonttype": 42}):
        fig.savefig(png_path, dpi=320, bbox_inches="tight", pad_inches=0.04)
        fig.savefig(png_path.with_suffix(".pdf"), bbox_inches="tight", pad_inches=0.04)
    plt.close(fig)
    return png_path


def _write_method(path, args, n_subjects):
    text = f"""# Intra-ROI vs Between-ROI FC Run-Baseline Comparison

This companion analysis uses the same ROI masks, voxel weights, and
intra-vs-between FC definitions as `med_effects.py`.

For each subject/session, FC was recomputed separately from beta run 1 and beta
run 2.

The plotted run-to-run baseline is a within-session variability magnitude:
`0.5 * (abs(OFF_run2 - OFF_run1) + abs(ON_run2 - ON_run1))`.

The plotted medication effect is recomputed as a matched run-level magnitude:
`0.5 * (abs(ON_run1 - OFF_run1) + abs(ON_run2 - OFF_run2))`.

These two subject-level quantities are compared with a paired one-sample test
on `run-level medication-effect magnitude - run-to-run baseline`, separately for
intra-ROI FC, between-ROI FC, and the primary intra-minus-between contrast.

Subjects included in the comparison: {n_subjects}.
Connectivity metric: `{M.INTRA_BETWEEN_FC_METRIC}`.
Voxel selection: `{args.voxel_selection}`.
"""
    path.write_text(text, encoding="utf-8")


def main():
    args = build_parser().parse_args()
    M.INTRA_BETWEEN_FC_METRIC = args.intra_between_fc_metric
    weight_img, rois, weighted_rois, roi_threshold = _prepare_rois(args)
    run_values, roi_values = _compute_run_level_values(args, weight_img, rois)
    run_deltas = _complete_run_deltas(run_values)
    comparison_values = _comparison_subject_values(run_values, run_deltas)
    if comparison_values.empty:
        raise ValueError("No complete subjects had OFF/ON run1/run2 values")

    source_stats = _source_stats(comparison_values)
    paired_stats = _paired_comparison_stats(comparison_values)
    summary = pd.concat([source_stats, paired_stats], ignore_index=True)
    cols = _metric_columns(run_values)

    out = args.out
    out.mkdir(parents=True, exist_ok=True)
    roi_values.to_csv(out / "intra_vs_between_fc_run_roi_values.csv", index=False)
    run_values.to_csv(out / "intra_vs_between_fc_run_session_values.csv", index=False)
    run_deltas.to_csv(out / "intra_vs_between_fc_run_subject_session_deltas.csv", index=False)
    comparison_values.to_csv(out / "intra_vs_between_fc_run_baseline_subject_values.csv", index=False)
    summary.to_csv(out / "intra_vs_between_fc_run_baseline_summary.csv", index=False)
    metadata = {
        "connectivity_metric": M.INTRA_BETWEEN_FC_METRIC,
        "voxel_selection": args.voxel_selection,
        "roi_percentile": float(args.roi_percentile),
        "weight_threshold": float(roi_threshold),
        "n_rois": int(len(rois)),
        "n_weighted_vigour_rois": int(len(weighted_rois)),
        "n_run_sessions": int(run_values.shape[0]),
        "n_run_session_deltas": int(run_deltas.shape[0]),
        "n_complete_subjects": int(comparison_values["subject"].nunique()),
        "run_baseline_definition": "0.5 * (abs(OFF_run2 - OFF_run1) + abs(ON_run2 - ON_run1))",
        "run_level_medication_effect_definition": "0.5 * (abs(ON_run1 - OFF_run1) + abs(ON_run2 - OFF_run2))",
        "primary_comparison": "run-level medication-effect magnitude minus run-to-run baseline",
    }
    (out / "intra_vs_between_fc_run_baseline_metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    _write_method(out / "intra_vs_between_fc_run_baseline_method.md", args, comparison_values["subject"].nunique())
    figure_path = _plot_comparison(comparison_values, summary, out, cols["plot_ylabel"])
    print(f"Saved {figure_path}")
    print(f"Saved {figure_path.with_suffix('.pdf')}")
    print(f"Saved {out / 'intra_vs_between_fc_run_session_values.csv'}")
    print(f"Saved {out / 'intra_vs_between_fc_run_subject_session_deltas.csv'}")
    print(f"Saved {out / 'intra_vs_between_fc_run_baseline_subject_values.csv'}")
    print(f"Saved {out / 'intra_vs_between_fc_run_baseline_summary.csv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
