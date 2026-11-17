# Where the code and the manuscript disagreed

The analysis scripts in the working folder had drifted from the submitted
text: parameters were changed while the paper was being written, some steps
were run interactively and never saved, and a few code paths still carried
settings from earlier drafts. This file records every discrepancy found while
assembling this repository, what the code did, what the paper says, and how it
was resolved.

The rule applied throughout: **the manuscript governs**, except where the
manuscript contradicts its own reported numbers. Those cases are flagged in
section B and need a decision from the authors.

---

## A. Code corrected to match the manuscript

These were resolved by changing the code. Nothing in the paper needs to move.

### A1. Bootstrap replicates: 200 in the code, 500 in the paper

`GBM更新.R` set `B <- 200` with a comment saying it could be lowered on slow
machines; Methods and Results both state 500 bootstrap replicates for the
optimism correction. The repository uses `N_BOOT_OPTIMISM <- 500`
(`R/utils.R`), and `R/11_internal_validation.R` applies it.

The same applies to the AUC comparison, which is 500 replicates throughout.

### A2. AUC comparison: DeLong in the code, bootstrap in the paper

`ROC Curves of Four IR Indices.R` compared indices with `pROC::roc.test(...,
method = "delong")`, and only for two of the three contrasts (SPISE vs TyG and
SPISE vs HOMA-BETA; the HOMA-IR comparison is missing from the script).

Methods state that model-based AUCs were compared **by bootstrap with 500
replicates**, Bonferroni-corrected across three comparisons. `R/07` implements
the bootstrap version and all three contrasts.

If the original DeLong p-values were what went into the manuscript, the two
approaches usually agree closely at this sample size, but the numbers should
be checked against the bootstrap output before submission.

### A3. GBM predictors: nine variables in the code, six in the paper

`GBM更新.R` defined

```r
predictors <- c("spise", "age", "gender", "bmi", "smoking",
                "eos", "lymph", "neut", "crp")
```

which is a nine-variable set including BMI and CRP. The paper's consensus set
is six variables: SPISE, lymphocyte, neutrophil, eosinophil, cancer history,
cardiovascular disease. The repository uses `features$consensus` from
`R/09_feature_selection.R`, and CRP is not read at all.

`crp` does not appear in any table or figure in the paper, so the nine-variable
version was almost certainly a leftover from an earlier run.

### A4. SMOTE commented out in the saved script

`5折cv.R` has

```r
#sampling = "smote",  # 处理类别不平衡（如哮喘率低）
```

commented out, yet the Results and the Figure 4C legend both describe
comparison-stage resampling. The live run clearly had it enabled — the
out-of-fold mean predicted probability is about 0.44, not 0.084, which is the
signature of a balanced training prior.

The repository enables `sampling = "smote"` in `R/10` and documents the
consequence in the README and in `R/15`.

### A5. Split described in a comment did not match the split used

`5折cv.R` carries the comment `# 时间拆分：2009-2014 训练，2015-2018 验证`, but
the objects passed to it (`data_train`, `data_extest`) were the 2009–2016 and
2017–2018 sets, which is what the paper reports. The comment is stale; the
split itself was correct. `R/02_build_cohort.R` now performs the split
explicitly by cycle and asserts the resulting 2445 / 2143.

### A6. Fold seeds did not match 5 × 10

`5折cv.R` passed `seeds = lapply(1:51, function(x) sample.int(100000, 300))`.
Fifty resamples need at least 51 seed slots, and tuning adds one more per
resample, so the extra tuning iterations fall back to R's global stream. The
fold structure was therefore reproducible across algorithms but not fully
specified.

`R/10` builds a seed vector of length `folds × repeats + 1`, with 200 integers
per resample, and reuses it for all seven algorithms so the comparison stays
paired.

### A7. Feature selection was never saved as a script

`Boruta()`, `cv.glmnet()` and `rfe()` do not appear in any `.R` file in the
working folder; the results survive only as objects inside `feature.rds`
(`boruta_feature`, `lasso_feature`, `rfe_features`, `consensus_feature`).
`R/09_feature_selection.R` reconstructs the three selections from the Methods
description.

The reconstruction is close but not identical, and it cannot be made identical
without the original seeds:

| Step | Published | Reconstructed (seed 42) |
|---|---|---|
| Boruta | 12 variables | **12, exactly the same set** |
| LASSO (lambda.min) | 11 variables | 10 — the same set minus neutrophil count |
| RF-RFE | optimum at 12 features | see `output/tables/Figure3D_rfe_profile.csv` |
| Consensus | 6 features | 5 — neutrophil count missing |

LASSO and RFE both draw a random fold split, and neither seed is recorded
anywhere, so the last variable in or out is luck. Neutrophil count sits exactly
on that boundary: over 40 different seeds the `lambda.min` solution never
retains it, while every other feature of the published set is stable.

Because the six-feature consensus drives every model in steps 10–14, the
repository archives it as `CONSENSUS_PUBLISHED` in `R/utils.R` and uses that
set downstream. `R/09` still runs all three algorithms, plots Figure 3 from the
fresh run, and prints how far the reconstruction is from the archive. This is
the only place in the repository where a result is taken from the archive
rather than recomputed, and it is flagged in the console output.

### A8. Hard-coded results inside plotting code

`rcs plot.R` hard-coded `geom_vline(xintercept = 4.55)` and
`p_nonlinear_txt <- paste0("<0.001")`, so the figure would keep showing the
published values even if the model changed. `R/05` takes the breakpoint from
the `segmented` fit and the p values from `anova()`.

### A9. A second, unrelated analysis inside `rcs plot.R`

The same file contains a full SPISE-and-osteoporosis analysis (IFG and diabetes
subgroups, `OST` as the outcome) belonging to a different project. It is not
part of this study and is not reproduced here.

### A10. Scripts depended on objects that only existed in the R session

Nearly every script assumed `data_all`, `data_train`, `data_extest`,
`final_gbm`, `nt`, `predictors` and friends were already in the global
environment, with no code producing them. The repository replaces this with
explicit `read_step()` / `write_step()` calls: each script reads the object the
previous one wrote and saves its own.

### A11. The candidate list had race in place of red cell count

Boruta's paragraph in the Results names the eighteen candidates implicitly, by
listing twelve retained and six rejected:

* retained — SPISE, age, fasting glucose, HbA1c, monocyte, lymphocyte,
  neutrophil, eosinophil, LDL cholesterol, total cholesterol, cancer, CVD
* rejected — alcohol, smoking, sex, red cell count, platelet count, RDW

Race is in neither list. The working script's candidate set had `RIDRETH1`
where red cell count (`LBXRBCSI`) should be, eighteen variables either way.
That single substitution changes the Boruta result: with race in, thirteen
variables are retained and five rejected, and Figure 3A no longer matches the
text. `R/09` uses the paper's list; race stays in the cohort because Models 2
and 3 adjust for it, but it is not offered to the selection step.

### A12. RF-RFE searched 2–18 features, the paper reports 4–16

`caret::rfe()` was called with `sizes = seq(2, length(candidates))` in the
reconstruction, i.e. 2 to 18. The Results state that subsets of 4 to 16
features were evaluated. `R/09` uses `sizes = 4:16`.

### A13. The algorithm that gets locked is the published one

`R/11` applies the prespecified rule — highest cross-validated AUC, with
algorithms equal at three decimals treated as tied and the Brier score
breaking the tie — exactly as the manuscript describes it. Applied to a fresh
run of step 10 it does not reproduce the published outcome: the top six
algorithms sit within one standard error of each other, and with this run's
fold splits random forest comes out ahead of GBM.

That matters more than one changed name. A random forest refitted on the whole
development set has an apparent AUC of 0.999; the model the paper describes has
an apparent AUC of 0.722 and a Brier score that matches the published external
validation to four decimals. Locking the forest would make every number in
steps 12 to 14 describe a different model from the one in the manuscript.

So, as with the consensus features (A7), the archived choice is used: the rule
is applied and its winner is written to
`output/tables/Table3_model_selection.csv`, but the model that is fitted,
bootstrapped, validated and explained is `ALGORITHM_PUBLISHED` (GBM) with the
tuning values step 10 found for it. The console output states both choices
whenever they differ.

---

## B. Manuscript statements that conflict with the data

These cannot be fixed in the code without breaking the reported numbers, so
they need an editorial decision.

**Status after the author's review (2026-11-18):**

* **B1, B3, B6 — resolved.** The manuscript (`asthma正文.docx`, backup
  `asthma正文_备份_B类修订前.docx`) was edited: the Methods now define asthma
  via MCQ010 *and* MCQ035 (current asthma), the Results say "current asthma"
  instead of "ever asthma", the SHAP sample is 200 throughout, and the Q4
  interval reads 0.29–0.85.
* **B2 — not a discrepancy.** The author confirmed that fasting insulin was
  part of the GLU file in 2009–2012 (`LBXIN` inside `GLU_F`/`GLU_G`) and only
  moved to its own INS component from 2013–2014. `R/01`/`R/02` were updated
  accordingly; rebuilding the cohort from the raw XPT files reproduces 4588
  participants with insulin for 4533 of them (98.8%).
* **B4, B5, B7, B8 — kept as published**, by the author's decision. The
  sections below describe what the code produces; the manuscript wording
  stands.

### B1. The asthma definition is current asthma, not lifetime diagnosis

Methods (Definitions of SPISE and asthma) say:

> Asthma was defined as a "yes" answer to "Has a doctor or other health
> professional ever told you that you have asthma?" (item MCQ010), which
> captures a lifetime diagnosis.

The analysis does not do that. In the analytic cohort, `asthma` is 1 only when
**both** `MCQ010 = 1` (ever diagnosed) **and** `MCQ035 = 1` (still has asthma).
The two definitions give different samples:

| Definition | Cases in the 4588-participant cohort |
|---|---|
| `MCQ010 = 1` only (lifetime) | 637 |
| `MCQ010 = 1` and `MCQ035 = 1` (current) | **413** |

Every reported number — 413 cases (9.0%), 205 in the development set, 208 in
the validation set — belongs to the current-asthma definition.

**Resolved in the manuscript:** the Methods sentence now reads

> Asthma was defined as a "yes" answer to "Has a doctor or other health
> professional ever told you that you have asthma?" (item MCQ010) together
> with a "yes" answer to "Do you still have asthma?" (item MCQ035), which
> captures current physician-diagnosed asthma.

and the Results phrase "(ever asthma)" now reads "(current asthma)".
`R/02_build_cohort.R` implements the current-asthma definition;
`asthma_from_items()` is the single place to change it.

### B2. The insulin subset — resolved, not a discrepancy

Results, on the index comparison:

> this analysis was restricted to the 4532 participants with available fasting
> insulin (98.8% of the full sample of 4588)

The initial concern was that `LBXIN` only appears in a dedicated INS file from
2013 onward. That is true but incomplete: in the 2009–2010 and 2011–2012
cycles fasting insulin is part of the glucose file (`GLU_F`, `GLU_G`), and the
author confirmed the 4532 figure comes from the full five-cycle merge.

The repository was updated accordingly:

* `R/01` still downloads 58 files — `GLU_F`/`GLU_G` carry `LBXIN` alongside
  `LBXGLU`; only `INS_H`, `INS_I`, `INS_J` are separate downloads.
* `R/02` reads `LBXIN` from the GLU file for 2009–2012 and from the INS file
  for 2013–2018, and stores it in the cohort.
* `R/07` no longer re-reads the INS files; it uses the cohort's `LBXIN`
  column directly.

Rebuilding the cohort from the raw XPT files reproduces 4588 participants, of
whom **4533 (98.8%)** have insulin — one more than the 4532 quoted, a
difference of one fasting-subsample record that does not change any reported
number. Applying the additional `LBXGLU > 63` filter that `HOMA-BETA` needs,
the index comparison runs on 4524 participants with 406 events
(`output/tables/index_comparison_sample.csv`).

### B3. SHAP sample size: 200 in one paragraph, 100 in another

Methods says SHAP was applied "in 200 randomly selected participants from the
external validation set". Results said "the SHAP estimates rest on 100
participants". The code used 100.

**Resolved (author decision: 200).** The Results sentence now reads 200, and
`SHAP_N_EXPLAIN` in `R/utils.R` is 200; Figure 5 was regenerated on the
200-participant sample.

### B4. Python is cited but not used

Methods: "Analyses used R 4.4.2 and Python 3.13." No Python appears anywhere in
the analysis; SHAP is computed with the R package `kernelshap`. Unless Python
was used for something not represented in these scripts, the clause can be
dropped.

### B5. The unweighted sensitivity analysis does not reproduce

Results, sensitivity analyses:

> Repeating the per-unit analysis without survey weights attenuated the
> estimates towards the null, but they remained statistically significant
> across all three models (OR 0.99, 95% CI 0.98–0.99, p < 0.001;
> Supplementary Table S2)

The same model, refitted from the same specification without weights, gives
OR 0.85 (95% CI 0.80–0.90) in Model 3, and 0.83 (0.79–0.87) and 0.82
(0.77–0.86) in Models 1 and 2 — see
`output/tables/TableS2_sensitivity_unweighted.csv`. Three things do not fit:

* **Magnitude.** 0.99 against 0.85. The weighted estimate is 0.89, so 0.99
  would mean the covariate-adjusted association almost disappears once the
  weights come off, which the data do not show.
* **Precision.** A 95% CI of 0.98–0.99 is 0.01 wide, which implies a standard
  error of about 0.003 on the log-odds scale. With 413 events the smallest
  attainable SE in this model is around 0.03 — ten times larger. That interval
  cannot come from this sample.
* **Direction.** The unweighted estimate is 0.85 against a weighted 0.89, i.e.
  slightly further from the null, not attenuated towards it.

An OR of 0.99 with a 0.01-wide interval looks like a variable on a much coarser
scale than one SPISE unit, so this is most likely a transcription from the
wrong row. Suggested replacement, if the code is right:

> Repeating the per-unit analysis without survey weights left the estimate
> essentially unchanged (OR 0.85, 95% CI 0.80–0.90, p < 0.001;
> Supplementary Table S2), which does not suggest a survey-design artefact.

### B6. The Q4 confidence interval in Table S1 is not symmetric

Results report the top quartile as OR 0.50 (95% CI 0.25–0.85, p = 0.014).
The point estimate and the p value reproduce exactly, but the interval does
not: on the log scale 0.50 sits 0.693 below the upper bound and 0.531 above the
lower one, and a Wald interval has to be symmetric. The repository computes
0.50 (0.29–0.85), which is symmetric.

The lower bound in the paper should read 0.29. **Resolved in the manuscript:**
the Results now read "OR 0.50, 95% CI 0.29–0.85, p = 0.014". Nothing else in
the sentence changes.

### B7. Table 1's prose reports means under the label "median"

`Table 1.docx` is a median (IQR) table, and `R/03` reproduces it: age 52
(37–67) overall, 52 (37–67) without asthma, 54 (41–67) with. The Results
paragraph, however, quotes means:

| Sentence in the Results | Paper | Mean (SD) in the data | Median in the data |
|---|---|---|---|
| "Median age was similar … (54 vs 52 years)" | 54 vs 52 | 53.7 vs 52.1 | 54 vs 52 |
| "higher body mass index (33.2 vs 29.4)" | 33.2 vs 29.4 | 33.20 vs 29.43 | 31.3 vs 28.2 |
| "a lower SPISE index (5.20 vs 6.03)" | 5.20 vs 6.03 | 5.20 vs 6.03 | 4.79 vs 5.64 |
| "eosinophil (0.25 vs 0.19)" | 0.25 vs 0.19 | 0.25 vs 0.19 | 0.20 vs 0.20 |
| "triglycerides (115 vs 111)" | 115 vs 111 | 115.22 vs 110.55 | 101 vs 96 |

Every figure in the paragraph is a mean, not a median. `R/03` now writes both
summaries (`Table1_baseline.csv` for the table as published,
`Table1_baseline_means.csv` for the numbers the text quotes), so the wording
can be fixed without recomputing anything: "Median age" should read "Mean
age".

Two smaller mismatches in the same paragraph, both in the mean version:
neutrophil count is 4.39 in the asthma group against the 4.29 quoted, and
monocyte count is higher, not lower, in the asthma group (0.56 vs 0.53,
against the "0.51 vs 0.53" in the text). HDL cholesterol p is 0.706 here and
0.071 in the text, which looks like a dropped zero.

### B8. The external calibration slope of 2.472 does not fit the rest of the paragraph

Results, temporal validation:

> Calibration drifted markedly: the calibration slope was 2.472 […] The Brier
> score was 0.0862, close to the null-model value of 0.0876.

The same locked model, applied to the same 2143 participants, gives a
calibration slope of 0.699, a Brier score of 0.0861 and the same null-model
value of 0.0876. The Brier scores agree to four decimals, so the model and the
predictions are the same; the slope is the outlier.

A slope of 2.472 means the log-odds of the prediction have to be stretched by
a factor of 2.5 to match the observed risk. With predictions centred on 8.5%
against an observed 9.7% — which is what the data show — that stretch would
move the fitted probabilities far from the base rate and push the Brier score
above the null value, not to within 0.0014 of it. The Brier score and the slope
cannot both be right.

`R/12` writes the slope it computes
(`output/tables/Table4_external_validation.csv`); if 0.70 is the true value the
sentence should read that, and "drifted markedly" should probably become
"drifted": a slope of 0.70 means the predictions are too extreme, which is
mild miscalibration rather than a collapse.

### B9. The index comparison: fresh run finds SPISE vs TyG significant too

The manuscript says SPISE's advantage over the two HOMA indices remained
significant after Bonferroni correction, implying the comparison against TyG
did not. With the corrected insulin merge the comparison now runs on 4524
participants (406 events) across all five cycles, and the fresh run
(`output/tables/Table2_index_AUC_comparisons.csv`) gives:

| Comparison | ΔAUC (95% CI) | p (Bonferroni ×3) |
|---|---|---|
| SPISE vs TyG | 0.011 (0.004 to 0.018) | < 0.001 |
| SPISE vs HOMA-IR | 0.010 (0.002 to 0.018) | 0.024 |
| SPISE vs HOMA-BETA | 0.011 (0.004 to 0.018) | 0.012 |

All three differences are significant at the corrected 0.05 level, with TyG
the *most* significant — the opposite pattern to the manuscript's wording.
The point estimates line up with the paper (SPISE highest, the other three
within 0.001 of each other), so the models are the same; the significance
pattern does not reproduce. The earlier code compared AUCs only for SPISE vs
TyG and SPISE vs HOMA-BETA (`KNOWN_DIFFERENCES` A2), and its seeds were not
recorded, so the exact earlier result cannot be recovered. Per the author's
decision the manuscript wording stands; the code reports what it computes.

---

## C. Scripts not carried into the repository

| File | Why |
|---|---|
| `模型比较2.R` | earlier version: 8:2 single split instead of repeated CV, and nine algorithms including RPART, SVM, NNET and XGBoost. Superseded by `5折cv.R` and not used for any reported result. |
| `baseline_asthma.R` | duplicate of `清洗+基线.R` with a renamed copy of the same description function. |
| `forestplot.R` / `forestplot_asthma.R` | identical files; the forest plot is regenerated in `R/06`. |
| `re_calibrated_curves.R`, `calibration curves.R`, `re_ROC_curves.R`, `re_dca_curves.R`, `roc_curves.R`, `dca_curves.R` | one-off figure variants that read objects from the session; their content is folded into `R/12` and `R/13`. |
| `Bootstrap可视.R`, `shap.R`, `shap组合.R`, `SHAP Dependence of SPISE.R`, `Boruta plot.R`, `Venn.R`, `bump_asthma_new.R`, `甜甜圈asthma.R`, `corplot.R`, `save_pred.R`, `AUC点估计.R`, `clinpubr_asthma.R`, `dca可用阈值.R` | figure grooming and intermediate checks, each rewritten as a section of the corresponding analysis script. |
| `spise.R`, `清洗+基线.R`, `代码.R` | the data assembly these implement is consolidated in `R/02`, `R/03` and `R/04`. |

The originals are kept in the working folder next to this repository; they are
not part of the public release because they cannot run without the session
state they were written in.

---

## D. Reporting conventions

Small things that were inconsistent in the draft and are now uniform in the
output tables:

* p values as `p = 0.047` or `p < 0.001`, italic p, space on both sides.
* Odds ratios as `0.89 (0.81–0.98)` with an en dash.
* Thousands separators from five digits up (BMC style): `4588` but `461,000`.
* `NA` shown as an empty cell in the CSV tables rather than the string "NA".
