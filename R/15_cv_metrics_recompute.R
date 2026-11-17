# =============================================================================
# 15_cv_metrics_recompute.R -- put the cross-validated metrics on the
#                             population scale
#
# SMOTE balances every training fold to 1:1, so the out-of-fold probabilities
# that step 10 stores live on a 50% prior. Their AUC is unaffected (a monotone
# transformation leaves the ranking alone) but their Brier score is not
# comparable with a value computed at the population event rate of about 8-9%.
#
# The shift is applied on the log-odds scale:
#     logit(p_population) = logit(p_balanced) + [logit(prev) - logit(0.5)]
# which changes the intercept only, not the slope or the AUC.
#
# Two summaries are produced:
#   - per repeat: metrics computed inside each of the 10 repetitions, then
#     averaged (mean +/- SD across repetitions)
#   - pooled: metrics computed once on the repeat-averaged probabilities
#
#   Rscript R/15_cv_metrics_recompute.R
# =============================================================================

source("R/utils.R")
need("dplyr", "pROC")

oof <- read_step("cv_predictions_oof")

prev <- mean(oof$y)
shift <- prior_shift(prev)

cat(sprintf("Event rate in the development set: %.4f\n", prev))
cat(sprintf("Null-model Brier: %.4f\n", prev * (1 - prev)))
cat(sprintf("Prior shift applied to the log-odds: %.3f\n", shift))

metric_one <- function(p, y) {
  p <- clip_prob(p)
  c(AUC = as.numeric(pROC::auc(pROC::roc(y, p, quiet = TRUE, direction = "<"))),
    Brier = brier(y, p),
    Slope = calibration_slope(y, p),
    Intercept = calibration_intercept(y, p),
    MeanP = mean(p))
}

# ---------------------------------------------------------------------------
# 1. Per repetition, then averaged
# ---------------------------------------------------------------------------
per_rep <- oof %>%
  dplyr::group_by(Model, Resample) %>%
  dplyr::group_modify(function(.x, .y) {
    m <- metric_one(.x$prob, .x$y)
    p_shifted <- stats::plogis(logit(.x$prob) + shift)
    data.frame(as.list(m),
               Brier_corrected = brier(.x$y, p_shifted),
               MeanP_corrected = mean(p_shifted))
  }) %>%
  dplyr::ungroup()

summ <- per_rep %>%
  dplyr::group_by(Model) %>%
  dplyr::summarise(
    n_repeat = dplyr::n(),
    dplyr::across(c(AUC, Brier, Brier_corrected, Slope, MeanP_corrected),
                  list(mean = mean, sd = stats::sd),
                  .names = "{.col}_{.fn}"),
    .groups = "drop"
  ) %>%
  dplyr::arrange(dplyr::desc(AUC_mean))

cat("\n--- per repetition, then averaged (mean +/- SD across ",
    CV_REPEATS, " repetitions) ---\n", sep = "")
print(as.data.frame(summ), row.names = FALSE)
save_table(summ, "Table3_cv_performance_corrected.csv")
save_table(per_rep, "TableS_cv_by_repeat_corrected.csv")

# ---------------------------------------------------------------------------
# 2. Pooled on repeat-averaged probabilities
# ---------------------------------------------------------------------------
pooled <- oof %>%
  dplyr::group_by(Model, rowIndex) %>%
  dplyr::summarise(prob = mean(prob), y = dplyr::first(y), .groups = "drop") %>%
  dplyr::group_by(Model) %>%
  dplyr::group_modify(function(.x, .y) {
    m <- metric_one(.x$prob, .x$y)
    p_shifted <- stats::plogis(logit(.x$prob) + shift)
    data.frame(as.list(m),
               Brier_corrected = brier(.x$y, p_shifted),
               MeanP_corrected = mean(p_shifted))
  }) %>%
  dplyr::ungroup() %>%
  dplyr::arrange(dplyr::desc(AUC))

cat("\n--- pooled over repetitions ---\n")
print(as.data.frame(pooled), row.names = FALSE)
save_table(pooled, "Table3_cv_performance_pooled.csv")

# ---------------------------------------------------------------------------
# 3. Reference values for the manuscript
# ---------------------------------------------------------------------------
cat("\n--------------------------------------------------------------------\n")
cat("When reading these numbers remember two things:\n")
cat("  1. Brier (raw) is on the balanced 50% scale and is meaningful only for\n")
cat("     ranking algorithms trained the same way.\n")
cat("  2. Brier_corrected is on the population scale and can be compared with\n")
cat("     the null-model value above. A Brier score below the null value means\n")
cat("     the model carries information; the gap is usually small at this\n")
cat("     event rate because most participants are true negatives.\n")
cat("  3. A low Brier score alone does not indicate good calibration: a model\n")
cat("     whose probabilities stay close to the base rate can score well while\n")
cat("     its calibration slope is far from 1. Read Brier together with Slope.\n")
cat("--------------------------------------------------------------------\n")

stamp_session("15_cv_metrics")
