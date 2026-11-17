# SPISE and asthma in US adults: analysis code

Analysis code for the study *"The Single Point Insulin Sensitivity Estimator
(SPISE) and asthma in US adults: a cross-sectional analysis and prediction
model with temporal external validation"*.

The pipeline takes the public NHANES files for 2009–2018, builds the analytic
cohort, and reproduces every number, table and figure in the manuscript:
survey-weighted association models, restricted cubic splines and the
breakpoint, subgroup analyses, the comparison of four insulin-sensitivity
indices, consensus feature selection, the seven-algorithm comparison under
repeated cross-validation, internal validation, temporal external validation
with recalibration, decision curve analysis and SHAP attribution.

---

## Data

All data are public and are not distributed with this repository.

| Source | Files | Notes |
|---|---|---|
| NHANES 2009–2010 to 2017–2018 | 58 `.XPT` files | `R/01_download_data.R` fetches them into `data/raw/` |

Download page: <https://wwwn.cdc.gov/nchs/nhanes/continuousnhanes/default.aspx>

If automated download fails (the CDC server throttles repeated requests),
download the files for each cycle by hand and place them in `data/raw/` with
their original names (`DEMO_F.XPT`, `CBC_H.XPT`, …). The script skips files
that are already present, so it can simply be re-run.

Components used: `DEMO`, `BMX`, `TRIGLY`, `HDL`, `TCHOL`, `GLU`, `GHB`, `CBC`,
`MCQ`, `SMQ`, `ALQ`, and `INS` (fasting insulin, available only from 2013
onward).

---

## Quick start

```bash
git clone <this-repo> && cd nhanes-spise-asthma

Rscript R/00_setup.R            # checks packages and prints the settings
Rscript R/01_download_data.R    # fetches the NHANES files (~120 MB)
Rscript R/02_build_cohort.R     # builds the analytic cohort

Rscript run_all.R --no-download # every remaining step, in order
```

A single step can be re-run on its own, because each one reads the object
written by the step before it:

```bash
Rscript R/10_model_comparison_cv.R
Rscript R/11_internal_validation.R
```

Runtimes on a 2020 laptop (8 cores): steps 01–08 take a few minutes in total;
step 09 needs about 15 minutes (Boruta runs 200 forests and RF-RFE refits the
forest 5-fold for every candidate subset size); step 10 is the slow one at
roughly 20–40 minutes because seven algorithms are tuned under 5-fold ×
10-repeat cross-validation; steps 11–15 take about ten minutes together.

---

## What each script does

| Script | Analysis | Manuscript |
|---|---|---|
| `00_setup.R` | package check, prints the analysis constants | — |
| `01_download_data.R` | downloads the NHANES `.XPT` files | — |
| `02_build_cohort.R` | merges components, derives SPISE / asthma / covariates, applies eligibility, splits by cycle | Methods, Figure 2 |
| `03_baseline_table.R` | baseline characteristics, Mann–Whitney and chi-square | Table 1 |
| `04_weighted_logistic.R` | survey-weighted models 1–3 | Table 2 |
| `05_rcs_breakpoint.R` | knot selection by AIC, spline fit, segmented regression | Figure 1A, Results |
| `06_subgroup.R` | seven interaction tests, stratum-specific ORs | Figure 1C |
| `07_index_comparison.R` | SPISE vs TyG, HOMA-IR, HOMA-BETA; bootstrap AUC comparison | Figure 1B, Table 2 |
| `08_sensitivity.R` | unweighted model, SPISE quartiles, all three adjustment sets | Supplementary Tables S1, S2 |
| `09_feature_selection.R` | Boruta, LASSO, RF-RFE, consensus | Figure 3 |
| `10_model_comparison_cv.R` | seven algorithms, repeated CV with SMOTE | Figure 4A–C, Table 3 |
| `11_internal_validation.R` | model selection, bootstrap optimism | Figure S1 |
| `12_external_validation.R` | temporal validation, cross-fitted recalibration | Figure 4D–F, Table 4 |
| `13_decision_curve.R` | net benefit before and after recalibration | Figure 4G |
| `14_shap.R` | kernel SHAP attribution | Figure 5, Table 5 |
| `15_cv_metrics_recompute.R` | puts the cross-validated metrics on the population scale | Table 3 footnote, Figure 4C legend |

---

## Two things worth knowing before reading the code

**1. SMOTE in the comparison stage only.** The seven algorithms are compared
with SMOTE applied inside each training fold, which is the fair way to compare
algorithms that react differently to a 1:8 class imbalance. The consequence is
that the out-of-fold probabilities sit on a balanced (50%) prior, so their
Brier score is not comparable with the population event rate of ~8.4%.

The locked model is therefore refitted on the whole development set *without*
resampling, which is why its predicted probabilities — and the external
validation — are on the population scale. `R/15_cv_metrics_recompute.R`
shifts the cross-validated probabilities back by
`logit(0.084) − logit(0.5) = −2.39` on the log-odds scale before any absolute
comparison is made. AUC is unaffected by this shift; the Brier score is.

**2. A low Brier score is not the same as good calibration.** At a 8.4% event
rate a model whose probabilities stay near the base rate scores well on Brier
while its calibration slope is far from 1. This is why the Brier score is used
only to separate algorithms with equal discrimination, and why the calibration
slope is reported alongside it. `R/15` prints the null-model Brier (0.0769) as
the reference point.

---

## Expected results

The cohort-level and association results reproduce exactly. Everything that
depends on a random split, a bootstrap, or a machine-learning fit lands close
to the published value but not on it; the size of that drift is listed below so
it is not mistaken for a bug.

| Quantity | Manuscript | Fresh run |
|---|---|---|
| Analytic cohort | 4588, 413 with asthma (9.0%) | same |
| Development / validation | 2445 (205) / 2143 (208) | same |
| SPISE breakpoint | 4.55 (3.58–5.51) | 4.55 (3.59–5.51) |
| Model 3 OR per SPISE unit | 0.89 (0.81–0.98) | 0.89 (0.81–0.98) |
| Q4 vs Q1, fully adjusted | 0.50 (0.25–0.85) | 0.50 (0.29–0.85) |
| Boruta retained | 12 | 12, same set |
| RF-RFE optimum | 12 features | 12 features |
| Consensus | 6 features | archived (see A7) |
| Cross-validated AUC, GBM | 0.673 | 0.677 |
| Cross-validated Brier, GBM | 0.139 | 0.181 |
| Apparent AUC | 0.722 | 0.745 |
| Optimism | 0.031 | 0.036 |
| Optimism-corrected AUC | 0.691 | 0.710 |
| External AUC | 0.655 (0.615–0.695) | 0.650 (0.610–0.690) |
| External Brier (null) | 0.0862 (0.0876) | 0.0862 (0.0876) |
| External calibration slope | 2.472 | 0.679 |
| Recalibrated slope / Brier / AUC | 0.953 / 0.0851 / 0.647 | 0.960 / 0.0856 / 0.643 |
| Net benefit at 10%, original | 0.016 | 0.020 |
| Net benefit at 10%, treat all | −0.003 | −0.003 |
| SHAP share of SPISE | 43.8% | 26.5% |

The pipeline is deterministic given `SEED`: two runs of the same code give the
same output. It is not bit-identical to the published analysis, because the
original fit was produced with seeds that were never recorded and with
selection results that are archived rather than recomputed (A7). The external
Brier score and its null value are the sharpest evidence that the locked model
is essentially the published one; the calibration slope and the SHAP shares are
the two quantities that move most.

`R/02_build_cohort.R` warns if the cohort size does not match, because a
mismatch there means the raw files or the eligibility rules differ from the
ones used in the paper.

---

## Environment

R 4.4.2 or later. Install the dependencies with:

```r
install.packages(c(
  "haven", "dplyr",
  "survey", "rms", "segmented", "pROC",
  "caret", "glmnet", "ranger", "gbm", "klaR", "e1071",
  "Boruta", "kernelshap", "shapviz",
  "ggplot2", "ggvenn",
  "doParallel", "foreach"
))
```

Decision curve analysis is implemented from its definition in `R/13` rather
than through a package, so the repository needs no dependency outside CRAN.

Every script writes `sessionInfo()` to `logs/`, so the exact package versions
used for a given run are recorded alongside the results.

---

## Reproducibility notes

* A single seed (`SEED` in `R/utils.R`) drives every stochastic step.
* The cross-validation fold structure is generated once and reused for all
  seven algorithms, so the comparison is paired.
* All analysis constants that also appear in the Methods — fold count,
  repeats, bootstrap replicates, Bonferroni denominators, SHAP sample size —
  are defined once at the top of `R/utils.R`. If a number in the manuscript
  changes, it changes there.
* Intermediate objects in `data/` are named after the step that produced them;
  deleting one forces only the dependent steps to re-run.
* Step 09 re-runs all three selection algorithms but uses the six-feature
  consensus archived in `R/utils.R` (`CONSENSUS_PUBLISHED`) for everything
  downstream. The original LASSO and RF-RFE fold seeds were never recorded and
  cannot be recovered, so the reconstruction reproduces Boruta exactly and the
  other two to within one variable. See `docs/KNOWN_DIFFERENCES.md`, A7.
* Step 11 applies the same rule to the algorithm choice and reports its
  winner, but fits and locks the published algorithm (`ALGORITHM_PUBLISHED`,
  GBM) whenever the two disagree — a random forest locked instead would have
  an apparent AUC of 0.999 and would not be the model the paper describes.
  See `docs/KNOWN_DIFFERENCES.md`, A13.

---

## Known differences between this code and the manuscript

`docs/KNOWN_DIFFERENCES.md` lists every place where the original analysis
scripts and the submitted text disagreed, and how each was resolved. Section A
lists thirteen code-side fixes. Section B lists eight statements in the manuscript
that the data do not support — including the asthma definition, the size of the
insulin subset, two numbers in the sensitivity analyses, and the external
calibration slope — and those need an editorial decision before submission.

---

## Citation

If you use this code, please cite the paper (details on acceptance) and the
NHANES data:

> National Center for Health Statistics. National Health and Nutrition
> Examination Survey, 2009–2018. Hyattsville, MD: NCHS.

---

## Licence

MIT. See `LICENSE`. The NHANES data carry their own terms of use.
