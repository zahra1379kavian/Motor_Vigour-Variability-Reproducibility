#!/usr/bin/env python3
"""Check Supplementary Figure 4 FEAT provenance files."""

from __future__ import annotations

from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
FIGURE = ROOT / "figures" / "supplementary" / "supp_figure_04_standard_glm_feat_brain_slices.png"
RESULTS = ROOT / "results" / "supplementary" / "figure_04_standard_glm_feat"
REQUIRED = (
    FIGURE,
    RESULTS / "feat_design" / "design.fsf",
    RESULTS / "feat_design" / "cope1_design.fsf",
    RESULTS / "feat_stats" / "thresh_zstat1_z3p1.nii.gz",
    RESULTS / "feat_stats" / "cluster_zstat1_z3p1_std.txt",
    RESULTS / "provenance" / "mixed_model.fsf",
    ROOT / "analysis" / "fsl_glm" / "run_mixed_from_template.py",
)


def main() -> int:
    missing = [path for path in REQUIRED if not path.exists()]
    if missing:
        print("Missing Supplementary Figure 4 provenance files:")
        for path in missing:
            print(f"- {path}")
        return 1
    print(f"Supplementary Figure 4 figure: {FIGURE}")
    print(f"Supplementary Figure 4 provenance: {RESULTS}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
