#!/usr/bin/env python3
"""Paper-ready first-row medication FC change-consistency figure.

This uses the existing summary CSVs from
figures/med_effects/med_fc_change_consistency and redraws only the
vigour-weighted Pearson and mutual-information panels.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.text import Text
import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
DEFAULT_INPUT_DIR = ROOT / "results" / "supplementary" / "figure_09_medication_fc_consistency"
DEFAULT_OUTPUT_STEM = DEFAULT_INPUT_DIR / "med_fc_change_consistency_first_row_pearson_mi"

FC_METRICS = ("pearson_r", "mutual_info_quantile")
X_POS = np.array([0.0, 0.58])


def _load_inputs(input_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    noise_path = input_dir / "change_vs_noise_subject.csv"
    summary_path = input_dir / "change_consistency_summary.csv"
    if not noise_path.exists() or not summary_path.exists():
        raise RuntimeError(f"Missing consistency CSVs under {input_dir}")
    return pd.read_csv(noise_path), pd.read_csv(summary_path)


def _draw_metric_axis(
    ax: plt.Axes,
    noise: pd.DataFrame,
    summary: pd.DataFrame,
    fc_metric: str,
    *,
    show_ylabel: bool,
) -> list:
    d = noise.loc[noise["fc_metric"] == fc_metric]
    if d.empty:
        raise RuntimeError(f"No subject rows found for {fc_metric}")
    rows = summary.loc[summary["fc_metric"] == fc_metric]
    if rows.empty:
        raise RuntimeError(f"No summary row found for {fc_metric}")
    g = rows.iloc[0]

    within = d["within_state_dist"].to_numpy()
    between = d["between_state_dist"].to_numpy()
    for wi, bi in zip(within, between):
        ax.plot(X_POS, [wi, bi], color="0.72", lw=1.15, zorder=1)

    within_points = ax.scatter(
        np.full_like(within, X_POS[0], dtype=float),
        within,
        color="#2c7fb8",
        s=38,
        zorder=2,
        label="within-state run 1 vs 2",
    )
    between_points = ax.scatter(
        np.full_like(between, X_POS[1], dtype=float),
        between,
        color="#f03b20",
        s=38,
        zorder=2,
        label="ON vs OFF",
    )

    ax.set_title(
        f"d={g['excess_cohens_d']:.2f}, p={g['excess_p_signflip']:.3g}",
        fontsize=14,
        pad=8,
    )
    ax.set_xticks(X_POS)
    ax.set_xticklabels(["run 1 vs run 2", "ON vs OFF"], fontsize=13)
    ax.set_xlim(X_POS[0] - 0.10, X_POS[1] + 0.10)
    ax.tick_params(axis="y", labelsize=12)
    ax.tick_params(axis="x", pad=7)
    if show_ylabel:
        ax.set_ylabel("Correlation distance", fontsize=14)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    return [within_points, between_points]


def make_figure(input_dir: Path, output_stem: Path) -> None:
    noise, summary = _load_inputs(input_dir)
    missing = [metric for metric in FC_METRICS if metric not in set(summary["fc_metric"])]
    if missing:
        raise RuntimeError(f"Missing requested metrics in summary CSV: {', '.join(missing)}")

    fig, axes = plt.subplots(1, len(FC_METRICS), figsize=(8.2, 4.2), squeeze=False)
    for idx, fc_metric in enumerate(FC_METRICS):
        _draw_metric_axis(
            axes[0, idx],
            noise,
            summary,
            fc_metric,
            show_ylabel=(idx == 0),
        )

    output_stem.parent.mkdir(parents=True, exist_ok=True)
    fig.subplots_adjust(left=0.10, right=0.99, bottom=0.12, top=0.91, wspace=0.16)
    for text in fig.findobj(match=Text):
        text.set_fontweight("bold")
    fig.savefig(output_stem.with_suffix(".png"), dpi=300, bbox_inches="tight", pad_inches=0.02)
    fig.savefig(output_stem.with_suffix(".pdf"), bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, default=DEFAULT_INPUT_DIR)
    parser.add_argument("--output-stem", type=Path, default=DEFAULT_OUTPUT_STEM)
    args = parser.parse_args()
    make_figure(args.input_dir, args.output_stem)
    print(f"Wrote {args.output_stem.with_suffix('.png')}")
    print(f"Wrote {args.output_stem.with_suffix('.pdf')}")


if __name__ == "__main__":
    main()
