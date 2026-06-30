#!/usr/bin/env python3
"""Regenerate supplementary Figure 16: task-activation GVS feature deltas."""

from __future__ import annotations

import shutil
from pathlib import Path

from _runner import run_module


ROOT = Path(__file__).resolve().parents[1]
RESULTS_DIR = ROOT / "results" / "supplementary" / "figure_16_gvs_task_activation_feature_deltas"
FIGURES_DIR = ROOT / "figures" / "supplementary"


def main() -> None:
    run_module(
        "plot_gvs_task_feature_delta_violins",
        ["plot_gvs_task_feature_delta_violins.py", "--in-dir", str(RESULTS_DIR)],
    )
    shutil.copy2(
        RESULTS_DIR / "gvs_vs_sham_feature_delta_violin_boxplot_session1_off.png",
        FIGURES_DIR / "supp_figure_16a_gvs_task_activation_feature_delta_session1_off.png",
    )
    shutil.copy2(
        RESULTS_DIR / "gvs_vs_sham_feature_delta_violin_boxplot_session2_on.png",
        FIGURES_DIR / "supp_figure_16b_gvs_task_activation_feature_delta_session2_on.png",
    )


if __name__ == "__main__":
    main()
