#!/usr/bin/env python3
"""Regenerate Figure 7A: GVS connectogram for vigour-network ROI betas."""

from _runner import run_script


if __name__ == "__main__":
    run_script("analysis/gvs_connectivity_coactivation/run_edge_connectivity_metric_sensitivity.py")
    run_script("analysis/gvs_connectivity_coactivation/plot_fdr_significant_edge_connectograms.py")
