# Code to manuscript map

Where each reported quantity is produced. Section names refer to the headings
in the manuscript; the script column is the file in `R/`.

## Methods-level quantities

| Quantity in the paper | Script | Output |
|---|---|---|
| Sample selection counts (Figure 2 flow) | `02_build_cohort.R` | printed to console, `data/cohort_analytic.csv` |
| Cohort size and event count | `02` | console check against 4588 / 413 |
| Development / validation split | `02` | `data/split.rds` |
| Knot number chosen by AIC | `05_rcs_breakpoint.R` | `output/tables/TableS_rcs_knot_selection.csv` |
| Breakpoint and its CI | `05` | `output/tables/TableS_breakpoint.csv` |
| Seven interaction p values | `06_subgroup.R` | `output/tables/TableS_interaction_tests.csv` |
| Bootstrap AUC contrasts | `07_index_comparison.R` | `output/tables/Table2_index_AUC_comparisons.csv` |
| Consensus feature set | `09_feature_selection.R` | `data/features.rds` (`consensus`, `consensus_reconstructed`), `Figure3E_venn.png` |
| Cross-validated AUC, SD, Brier, slope | `10_model_comparison_cv.R` | `Table3_cv_performance.csv` |
| Optimism and corrected AUC | `11_internal_validation.R` | `Table3_internal_validation.csv` |
| Recalibration parameters | `12_external_validation.R` | `TableS_recalibration_parameters.csv` |
| Net benefit at 10% and 15% | `13_decision_curve.R` | `Table4_net_benefit.csv` |
| Threshold window beating both defaults | `13` | `Table4_dca_window.csv`, `Table4_dca_curve.csv` |
| SHAP share per feature | `14_shap.R` | `Table5_shap_contribution.csv` |
| Population-scale Brier after SMOTE | `15_cv_metrics_recompute.R` | `Table3_cv_performance_corrected.csv` |

## Tables

| Table | Script | File |
|---|---|---|
| Table 1 — baseline characteristics | `03` | `output/tables/Table1_baseline.csv` |
| Table 2 — ORs for SPISE (models 1–3) | `04` | `Table2_SPISE_odds_ratios.csv`, `Table2_full_model3.csv` |
| Table 2 — AUC of the four indices | `07` | `Table2_index_AUCs.csv` |
| Table 3 — algorithm comparison | `10`, `15` | `Table3_cv_performance*.csv` |
| Table 4 — external validation | `12` | `Table4_external_validation.csv` |
| Table 5 — SHAP attribution | `14` | `Table5_shap_contribution.csv` |
| Supplementary tables | all | `TableS_*.csv`, `TableS1_sensitivity_quartiles.csv`, `TableS2_sensitivity_unweighted.csv` |

## Figures

| Figure | Panel | Script | File |
|---|---|---|---|
| 1 | A — RCS curve | `05` | `Figure1A_rcs.png` |
| 1 | B — ROC of four indices | `07` | `Figure1B_roc_indices.png` |
| 1 | C — subgroup forest plot | `06` | `Figure1C_subgroups.png` |
| 2 | — participant flow | not generated in code | counts printed by `02` |
| 3 | A — Boruta importance | `09` | `Figure3A_boruta.png` |
| 3 | B — LASSO deviance | `09` | `Figure3B_lasso_deviance.png` |
| 3 | D — RF-RFE profile | `09` | `Figure3D_rfe.png` |
| 3 | E — Venn | `09` | `Figure3E_venn.png` |
| 3 | F — correlation | `09` | `Figure3F_correlation.csv` |
| 4 | A–C — model comparison | `10` | `Table3_cv_performance.csv` (plot in manuscript) |
| 4 | D — calibration, external | `12` | `Figure4D_calibration.png` |
| 4 | F — ROC, external | `12` | `Figure4F_roc_external.png` |
| 4 | G — decision curve | `13` | `Figure4G_dca.png` |
| 5 | A, B — SHAP beeswarm, bar | `14` | `Figure5A_shap_beeswarm.png`, `Figure5B_shap_bar.png` |
| 5 | C — dependence | `14` | `Figure5C_shap_dependence_spise.png` |
| 5 | D, E — individual explanations | `14` | `Figure5D_waterfall_*.png` |
| S1 | — optimism distribution | `11` | `FigureS1_optimism.png` |

Figure 2 is a hand-assembled flow diagram; the counts it displays come from
the console output of `02_build_cohort.R`, which prints the number retained at
each eligibility step.

## Numbers that appear in the Abstract

| Statement | Source |
|---|---|
| 4588 participants, 413 with asthma (9.0%) | `02` |
| SPISE associated with asthma, OR from model 3 | `04` |
| Breakpoint 4.55 (3.58–5.51) | `05` |
| AUC 0.711 for SPISE in the index comparison | `07` |
| External AUC 0.655 (0.615–0.695) | `12` |
| Calibration slope 2.472 → 0.953 | `12` |
| SPISE accounts for 43.8% of SHAP attribution | `14` |
