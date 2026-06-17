#!/usr/bin/env python3
"""Regenerate supplementary Figures 6 and 7: threshold robustness summaries."""

from _runner import run_module


if __name__ == "__main__":
    raise SystemExit(
        run_module(
            "threshold_robustness_voxel_network",
            [
                "threshold_robustness_voxel_network.py",
                "--out-base",
                "results/supplementary/figure_06_07_threshold_robustness/vigour_network_threshold_robustness",
            ],
        )
    )
