#!/usr/bin/env python3
"""Regenerate supplementary Figure 13: original GVS time courses."""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path

from _runner import run_module


ROOT = Path(__file__).resolve().parents[1]
RESULTS_DIR = ROOT / "results" / "supplementary" / "figure_13_gvs_original_signals"
FIGURES_DIR = ROOT / "figures" / "supplementary"
SOURCE_NAME = "original_signal_gvs_01_to_09_colored_subplots_vertically_offset_session1_off_session2_on"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--panel",
        choices=("all", "task_activation", "vigour_network"),
        default="all",
        help="Panel to regenerate. Defaults to both Supplementary Figure 13 panels.",
    )
    return parser.parse_args()


def _run_panel(panel_dir: str, figure_stem: str, extra_args: list[str]) -> None:
    out_dir = RESULTS_DIR / panel_dir
    run_module(
        "plot_gvs_original_signal_colored_subplots",
        [
            "plot_gvs_original_signal_colored_subplots.py",
            "--out-dir",
            str(out_dir),
            *extra_args,
        ],
    )
    shutil.copy2(out_dir / f"{SOURCE_NAME}.png", FIGURES_DIR / f"{figure_stem}.png")
    pdf_path = out_dir / f"{SOURCE_NAME}.pdf"
    if pdf_path.exists():
        shutil.copy2(pdf_path, FIGURES_DIR / f"{figure_stem}.pdf")


def main() -> None:
    args = parse_args()
    if args.panel in {"all", "task_activation"}:
        _run_panel(
            "task_activation",
            "supp_figure_13a_gvs_task_activation_original_signal",
            ["--task-map-bold-trials", str(ROOT / "data" / "external" / "Task_map_BOLD_trials.npy")],
        )
    if args.panel in {"all", "vigour_network"}:
        _run_panel("vigour_network", "supp_figure_13b_gvs_vigour_network_original_signal", [])


if __name__ == "__main__":
    main()
