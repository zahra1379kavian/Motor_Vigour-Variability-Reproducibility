#!/usr/bin/env python3
"""Regenerate supplementary Figure 17: voxel trial variability comparison."""

from __future__ import annotations

import shutil
from pathlib import Path

from _runner import run_module


ROOT = Path(__file__).resolve().parents[1]
RESULTS_DIR = ROOT / "results" / "supplementary" / "figure_17_voxel_trial_variability"
FIGURES_DIR = ROOT / "figures" / "supplementary"
SOURCE_STEM = "voxel_trial_variability_comparison_blue_orange_purple"
FIGURE_STEM = "supp_figure_17_voxel_trial_variability_comparison"


def main() -> None:
    run_module("plot_bold_voxel_trial_variability")
    shutil.copy2(RESULTS_DIR / f"{SOURCE_STEM}.png", FIGURES_DIR / f"{FIGURE_STEM}.png")
    shutil.copy2(RESULTS_DIR / f"{SOURCE_STEM}.pdf", FIGURES_DIR / f"{FIGURE_STEM}.pdf")


if __name__ == "__main__":
    main()
