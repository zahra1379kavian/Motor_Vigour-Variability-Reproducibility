#!/usr/bin/env python3
"""Regenerate Figure 7B: GVS connectogram for task-activation ROI betas."""

from _runner import run_script


if __name__ == "__main__":
    run_script("analysis/gvs_connectivity_coactivation/run_task_activation_connectogram_analysis.py")
