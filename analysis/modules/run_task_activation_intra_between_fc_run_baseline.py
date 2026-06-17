#!/usr/bin/env python3
"""Run-level intra-vs-between FC baseline for task-activation ROIs."""

from __future__ import annotations

import argparse
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
for path in (ROOT, HERE):
    if str(path) not in sys.path:
        sys.path.insert(0, str(path))

import med_effects as M  # noqa: E402
import med_effects_task_activation as task  # noqa: E402

RUN_BASELINE_SCRIPT = ROOT / "analysis" / "modules" / "run_intra_between_fc_run_baseline.py"
DEFAULT_OUT = ROOT / "results" / "supplementary" / "figure_12_run_baseline_fc" / "task_activation"


def _load_run_baseline_helpers():
    spec = importlib.util.spec_from_file_location("vigour_run_baseline", RUN_BASELINE_SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Could not import {RUN_BASELINE_SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


base = _load_run_baseline_helpers()


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--task-activation-map", type=Path, default=task.DEFAULT_TASK_ACTIVATION_MAP)
    parser.add_argument("--task-z-threshold", type=float, default=task.DEFAULT_TASK_Z_THRESHOLD)
    parser.add_argument("--session-manifest", type=Path, default=M.DEFAULT_SESSION_MANIFEST)
    parser.add_argument("--beta-root", type=Path, default=M.DEFAULT_BETA_ROOT)
    parser.add_argument("--subjects", nargs="+", default=None)
    parser.add_argument("--complete-subjects-only", action="store_true")
    parser.add_argument("--out", type=Path, default=DEFAULT_OUT)
    hemisphere_group = parser.add_mutually_exclusive_group()
    hemisphere_group.add_argument("--split-hemispheres", dest="split_hemispheres", action="store_true", default=True)
    hemisphere_group.add_argument("--no-split-hemispheres", dest="split_hemispheres", action="store_false")
    parser.add_argument("--exclude-rois", nargs="*", default=())
    parser.add_argument("--min-report-voxels", type=int, default=M.DEFAULT_MIN_REPORT_VOXELS)
    parser.add_argument("--min-lateralized-voxels", type=int, default=task.DEFAULT_MIN_LATERALIZED_VOXELS)
    parser.add_argument("--aal-version", default=M.DEFAULT_AAL_VERSION)
    parser.add_argument("--atlas-cache-dir", type=Path, default=M.DEFAULT_ATLAS_CACHE_DIR)
    parser.add_argument("--intra-between-fc-metric", choices=M.INTRA_BETWEEN_FC_METRICS, default=M.INTRA_BETWEEN_FC_METRIC)
    parser.add_argument("--mi-quantile-bins", type=int, default=M.DEFAULT_MI_QUANTILE_BINS)
    parser.add_argument("--random-state", type=int, default=0)
    return parser


def _missing_inputs(args):
    missing = task._missing_inputs(args)
    if not RUN_BASELINE_SCRIPT.exists():
        missing.append(f"{RUN_BASELINE_SCRIPT} (run-baseline helper script)")
    return missing


def _write_method(path: Path, args, n_subjects: int) -> None:
    text = f"""# Task-Activation Intra-ROI vs Between-ROI FC Run-Baseline Comparison

This companion analysis uses the task-activation ROI masks from
`med_effects_task_activation.py`: standard-GLM z-map voxels thresholded at
z >= {float(args.task_z_threshold):g}, unit voxel weights, and the same
lateralized AAL grouping used for the task-activation medication-change figure.

For each subject/session/run, intra-ROI voxel-pair FC and between-ROI FC were
computed separately.

The plotted run-to-run baseline is a within-session variability magnitude:
`0.5 * (abs(OFF_run2 - OFF_run1) + abs(ON_run2 - ON_run1))`.

The plotted medication effect is a matched run-level magnitude:
`0.5 * (abs(ON_run1 - OFF_run1) + abs(ON_run2 - OFF_run2))`.

The paired contrast tests `run-level medication-effect magnitude -
run-to-run baseline` across subjects, separately for intra-ROI FC,
between-ROI FC, and the intra-minus-between contrast.

Subjects included in the comparison: {n_subjects}.
Connectivity metric: `{M.INTRA_BETWEEN_FC_METRIC}`.
"""
    path.write_text(text, encoding="utf-8")


def main() -> int:
    args = build_parser().parse_args()
    M.INTRA_BETWEEN_FC_METRIC = args.intra_between_fc_metric
    base.M.INTRA_BETWEEN_FC_METRIC = args.intra_between_fc_metric
    missing = _missing_inputs(args)
    if missing:
        print("Missing required inputs:")
        for item in missing:
            print(f"- {item}")
        return 1

    task_img, task_values, _, rois, _ = task._prepare_task_rois(args)
    selected_task = np.isfinite(task_values) & (task_values >= float(args.task_z_threshold))

    run_values, roi_values = base._compute_run_level_values(args, task_img, rois)
    run_deltas = base._complete_run_deltas(run_values)
    comparison_values = base._comparison_subject_values(run_values, run_deltas)
    if comparison_values.empty:
        raise RuntimeError("No complete subjects had OFF/ON run1/run2 values")

    source_stats = base._source_stats(comparison_values)
    paired_stats = base._paired_comparison_stats(comparison_values)
    summary = pd.concat([source_stats, paired_stats], ignore_index=True)
    cols = base._metric_columns(run_values)

    out = args.out
    out.mkdir(parents=True, exist_ok=True)
    task._roi_summary(rois, args.task_z_threshold).to_csv(out / "task_activation_roi_definition.csv", index=False)
    roi_values.to_csv(out / "intra_vs_between_fc_run_roi_values.csv", index=False)
    run_values.to_csv(out / "intra_vs_between_fc_run_session_values.csv", index=False)
    run_deltas.to_csv(out / "intra_vs_between_fc_run_subject_session_deltas.csv", index=False)
    comparison_values.to_csv(out / "intra_vs_between_fc_run_baseline_subject_values.csv", index=False)
    summary.to_csv(out / "intra_vs_between_fc_run_baseline_summary.csv", index=False)
    metadata = {
        "connectivity_metric": M.INTRA_BETWEEN_FC_METRIC,
        "voxel_selection": task.TASK_VOXEL_SELECTION,
        "voxel_weighting": "unit_weights",
        "task_activation_map": str(args.task_activation_map),
        "task_activation_threshold": f"z >= {float(args.task_z_threshold):g}",
        "task_activation_voxels": int(np.count_nonzero(selected_task)),
        "n_rois": int(len(rois)),
        "n_run_sessions": int(run_values.shape[0]),
        "n_run_session_deltas": int(run_deltas.shape[0]),
        "n_complete_subjects": int(comparison_values["subject"].nunique()),
        "run_baseline_definition": "0.5 * (abs(OFF_run2 - OFF_run1) + abs(ON_run2 - ON_run1))",
        "run_level_medication_effect_definition": "0.5 * (abs(ON_run1 - OFF_run1) + abs(ON_run2 - OFF_run2))",
        "primary_comparison": "run-level medication-effect magnitude minus run-to-run baseline",
    }
    (out / "intra_vs_between_fc_run_baseline_metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    _write_method(out / "intra_vs_between_fc_run_baseline_method.md", args, comparison_values["subject"].nunique())
    figure_path = base._plot_comparison(comparison_values, summary, out, cols["plot_ylabel"])
    print(f"Saved {figure_path}")
    print(f"Saved {figure_path.with_suffix('.pdf')}")
    print(f"Saved {out / 'intra_vs_between_fc_run_session_values.csv'}")
    print(f"Saved {out / 'intra_vs_between_fc_run_subject_session_deltas.csv'}")
    print(f"Saved {out / 'intra_vs_between_fc_run_baseline_subject_values.csv'}")
    print(f"Saved {out / 'intra_vs_between_fc_run_baseline_summary.csv'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
