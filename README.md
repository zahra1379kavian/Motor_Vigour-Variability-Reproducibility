# Motor Vigour Variability Reproducibility

This repository contains the curated code, derived inputs, result tables, and
main figure exports for the motor-vigour variability paper. It is intentionally
smaller than the original working thesis repository: exploratory analyses and
large intermediate folders were left out.

Supplementary figures have not been added yet.

## Repository Layout

- `figures/main/`: publication-ready main figures as PNG and PDF.
- `results/main/`: machine-readable tables, JSON/NPZ summaries, and source
  figure exports grouped by paper figure.
- `data/`: compact derived maps, ablation inputs, metadata, and processed GVS
  connectivity tables needed by the included analyses.
- `analysis/modules/`: copied analysis modules from the working repository,
  with defaults updated for this curated layout.
- `analysis/gvs_connectivity_coactivation/`: GVS connectivity analysis code.
- `scripts/`: reviewer-facing entry points with descriptive figure names.

## Main Figures

| Paper panel | Curated figure file | Primary script |
| --- | --- | --- |
| Figure 2A | `figures/main/figure_02a_behavior_projection_subject_panel.png` | `scripts/figure_02a_behavior_projection.py` |
| Figure 2B | `figures/main/figure_02b_trial_variability_hypothesis.png` | `scripts/figure_02b_trial_variability_hypothesis.py` |
| Figure 3 | `figures/main/figure_03_vigour_network_voxel_weights.png` | `scripts/figure_03_vigour_network_map.py` |
| Figure 4 | `figures/main/figure_04_full_model_vs_task_only_anatomy.png` | `scripts/figure_04_ablation_anatomy.py` |
| Figure 5 | `figures/main/figure_05_full_model_vs_task_only_roi_summary.png` | `scripts/figure_05_ablation_roi_summary.py` |
| Figure 6A | `figures/main/figure_06a_medication_fc_vigour_network.png` | `scripts/figure_06a_medication_vigour_network_fc.py` |
| Figure 6B | `figures/main/figure_06b_medication_fc_task_activation.png` | `scripts/figure_06b_medication_task_activation_fc.py` |
| Figure 7A | `figures/main/figure_07a_gvs_connectogram_vigour_network.png` | `scripts/figure_07a_gvs_vigour_connectogram.py` |
| Figure 7B | `figures/main/figure_07b_gvs_connectogram_task_activation.png` | `scripts/figure_07b_gvs_task_activation_connectogram.py` |

See `figures/main/figure_manifest.csv` for the source result directory for
each panel.

## Python Environment

Create a local virtual environment and install the listed dependencies:

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

The code was validated with the thesis repository virtual environment.

## Reproducing Figures

Run scripts from the repository root:

```bash
python scripts/figure_03_vigour_network_map.py
python scripts/figure_04_ablation_anatomy.py
python scripts/figure_05_ablation_roi_summary.py
python scripts/figure_07a_gvs_vigour_connectogram.py
```

Figures 2A, 2B, 6A, 6B, and parts of 7B require subject-level beta,
behaviour, or atlas resources that are not packaged here because they are
large and may be access-controlled. The current figure exports and companion
result tables are included under `figures/main/` and `results/main/`.

## External Data

External data requirements are documented in `data/external/README.md`. The
packaged `data/` folder contains compact derived maps and processed tables, not
the full raw/preprocessed subject-level dataset.

