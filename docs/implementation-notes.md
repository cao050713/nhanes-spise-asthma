# Implementation notes

Decisions a reader needs in order to follow the pipeline alongside the
manuscript: which quantities are recomputed, which are archived, and how the
numbers in `output/` are formatted.

---

## 1. Archived constants

Two inputs of the published analysis depend on random fold splits whose seeds
were never recorded, so they cannot be regenerated bit for bit. Both are kept in
`R/utils.R` as named constants so that everything downstream of feature
selection is built on the same inputs as the reported analyses.

### `CONSENSUS_PUBLISHED` — the six-feature consensus

`R/09_feature_selection.R` runs all three selection algorithms from scratch and
writes their output to `output/tables/`
(`Figure3D_rfe_profile.csv`, the Boruta and LASSO summaries, and
`data/features.rds`, which stores both the fresh reconstruction and the archived
set). It also prints how far the reconstruction is from the archived six:
Boruta reproduces exactly, LASSO and RF-RFE differ by at most the last variable
in or out, because neither original seed is recoverable. `R/10`–`R/14` use
`CONSENSUS_PUBLISHED`.

### `ALGORITHM_PUBLISHED` — the locked learner

`R/11_internal_validation.R` applies the prespecified rule — highest
cross-validated AUC, algorithms equal at three decimals treated as tied and the
Brier score breaking the tie — and writes the winner to
`output/tables/Table3_model_selection.csv`. The model that is then bootstrapped,
validated and explained is `ALGORITHM_PUBLISHED` (GBM). The two are not
interchangeable: in a fresh run the top six algorithms sit within one standard
error of each other, while a random forest refitted on the whole development set
has an apparent AUC of 0.999 against the 0.722 of the locked model. The console
notes it whenever the rule and the archive pick different learners.

Everything else — cohort assembly, survey weights, splines and the breakpoint,
subgroup interaction tests, the index comparison, repeated cross-validation,
optimism correction, recalibration, decision curves and SHAP — is recomputed
from the raw NHANES files.

---

## 2. Choices that follow the manuscript rather than the defaults

* **Quartile cut points** in `R/08_sensitivity.R` come from the unweighted
  sample distribution of SPISE, which is what Table S1 reports; weighted
  quantiles (`svyquantile`) are not used anywhere.
* **Bootstrap replicates**: 500 throughout, for the optimism correction and for
  the AUC contrasts (Bonferroni-corrected across three comparisons).
* **SHAP sample**: 200 participants drawn from the external validation set
  (`SHAP_N_EXPLAIN`).
* **Baseline table**: `R/03_baseline_table.R` writes both the median (IQR)
  version that Table 1 reports and a mean (SD) version
  (`Table1_baseline_means.csv`), so that the numbers quoted in the Results can
  be checked against either summary.
* **Asthma definition**: current asthma, i.e. `MCQ010 = 1` and `MCQ035 = 1`.
  `asthma_from_items()` in `R/02_build_cohort.R` is the single place to change
  it; the derived variable is documented in `docs/data-dictionary.md`.

---

## 3. Reporting conventions in `output/`

* p values as `p = 0.047` or `p < 0.001`, italic p, one space on both sides.
* Odds ratios and AUC differences as `0.89 (0.81–0.98)`, with an en dash.
* Thousands separators from five digits up (BMC style): `4588` but `461,000`.
* Missing values are empty cells in the CSV tables, not the string `NA`.
