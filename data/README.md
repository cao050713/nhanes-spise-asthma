# data/

Nothing in this folder is committed to the repository.

| Path | Produced by | Contents |
|---|---|---|
| `raw/` | `R/01_download_data.R` | the NHANES `.XPT` files, one per component and cycle |
| `cohort.rds` | `R/02_build_cohort.R` | analytic cohort, one row per participant |
| `split.rds` | `R/02` | development and validation sets |
| `features.rds` | `R/09_feature_selection.R` | Boruta, LASSO, RF-RFE and consensus selections |
| `trained_models.rds` | `R/10_model_comparison_cv.R` | seven fitted `caret` objects (large, ~100 MB) |
| `cv_predictions_oof.rds` | `R/10` | out-of-fold predictions for every algorithm and repeat |
| `locked_model.rds` | `R/11_internal_validation.R` | the final GBM, its imputation parameters and tuning values |
| `external_predictions.rds` | `R/12_external_validation.R` | validation-set predictions, original and recalibrated |
| `shap_values.rds` | `R/14_shap.R` | kernel SHAP values |
| `cohort_analytic.csv` | `R/02` | plain-text copy of the cohort, for inspection |
| `file_manifest.csv` | `R/01` | which files were downloaded from where |

Deleting an object forces only the steps downstream of it to re-run. For
example, removing `trained_models.rds` costs another 20–40 minutes in `R/10`
but nothing upstream of it.
