#!/usr/bin/env python3
"""Regenerate supplementary Figure 10: Pearson GVS connectograms."""

from _runner import run_script


if __name__ == "__main__":
    run_script("analysis/gvs_connectivity_coactivation/run_edge_connectivity_metric_sensitivity.py")
    run_script("analysis/gvs_connectivity_coactivation/plot_fdr_significant_edge_connectograms.py")
    run_script(
        "analysis/gvs_connectivity_coactivation/run_task_activation_connectogram_analysis.py",
        [
            "run_task_activation_connectogram_analysis.py",
            "--metrics",
            "pearson_r",
            "mutual_info_quantile",
            "spearman_rho",
            "--reuse-trial-table",
        ],
    )
