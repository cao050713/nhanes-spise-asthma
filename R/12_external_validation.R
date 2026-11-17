# =============================================================================
# 12_external_validation.R -- temporal validation and recalibration (Figure 4D-F)
#
# The locked model is applied to 2017-2018 without retraining. Updating
# parameters (intercept and slope) are estimated by 5-fold cross-fitting
# inside the validation set, so that no participant contributes to both the
# estimation and the evaluation of the update.
#
#   Rscript R/12_external_validation.R
# =============================================================================

source("R/utils.R")
# caret supplies predict.train for the locked model object.
need("dplyr", "caret", "pROC", "ggplot2")

split  <- read_step("split")
locked <- read_step("locked_model")

val <- split$validation
y   <- as.integer(val$asthma)

# ---------------------------------------------------------------------------
# 1. Predictions from the locked model
# ---------------------------------------------------------------------------
# The locked fit is a caret train object (see R/11); predict() with
# type = "prob" applies the preprocessing stored inside it.
val_pred <- function(newdata) {
  as.numeric(stats::predict(locked$fit, newdata = newdata, type = "prob")[, "Yes"])
}
pred <- val_pred(val[, locked$vars])

stopifnot(length(pred) == length(y))

# ---------------------------------------------------------------------------
# 2. Performance before recalibration
# ---------------------------------------------------------------------------
roc_raw <- pROC::roc(y, pred, quiet = TRUE, direction = "<")
ci_raw  <- as.numeric(pROC::ci.auc(roc_raw))

orig <- data.frame(
  n = length(y), events = sum(y),
  AUC = as.numeric(pROC::auc(roc_raw)),
  AUC_lower = ci_raw[1], AUC_upper = ci_raw[3],
  slope = calibration_slope(y, pred),
  intercept = calibration_intercept(y, pred),
  Brier = brier(y, pred),
  Brier_null = brier_null(y),
  mean_predicted = mean(pred),
  observed = mean(y)
)
cat("\n--- locked model, external validation ---\n")
print(orig, row.names = FALSE)

# ---------------------------------------------------------------------------
# 3. Cross-fitted logistic recalibration
# ---------------------------------------------------------------------------
set.seed(SEED)
folds <- caret::createFolds(y, k = 5, list = TRUE)
lp <- logit(pred)
pred_upd <- numeric(length(y))
upd_params <- matrix(NA_real_, nrow = length(folds), ncol = 2,
                     dimnames = list(NULL, c("intercept", "slope")))

for (i in seq_along(folds)) {
  test  <- folds[[i]]
  train <- setdiff(seq_along(y), test)
  fit_upd <- stats::glm(y[train] ~ lp[train], family = binomial())
  upd_params[i, ] <- stats::coef(fit_upd)
  pred_upd[test] <- stats::plogis(stats::coef(fit_upd)[1] +
                                    stats::coef(fit_upd)[2] * lp[test])
}

cat("\nPer-fold updating parameters:\n"); print(round(upd_params, 3))

roc_upd <- pROC::roc(y, pred_upd, quiet = TRUE, direction = "<")
ci_upd  <- as.numeric(pROC::ci.auc(roc_upd))

upd <- data.frame(
  AUC = as.numeric(pROC::auc(roc_upd)),
  AUC_lower = ci_upd[1], AUC_upper = ci_upd[3],
  slope = calibration_slope(y, pred_upd),
  intercept = calibration_intercept(y, pred_upd),
  Brier = brier(y, pred_upd),
  mean_predicted = mean(pred_upd), observed = mean(y)
)
cat("\n--- after cross-fitted recalibration ---\n")
print(upd, row.names = FALSE)

comparison <- rbind(
  cbind(stage = "Original", orig[, c("AUC", "AUC_lower", "AUC_upper", "slope",
                                     "intercept", "Brier")]),
  cbind(stage = "Recalibrated", upd[, c("AUC", "AUC_lower", "AUC_upper", "slope",
                                        "intercept", "Brier")])
)
comparison$AUC_95CI <- sprintf("%.3f (%.3f-%.3f)", comparison$AUC,
                               comparison$AUC_lower, comparison$AUC_upper)
print(comparison[, c("stage", "AUC_95CI", "slope", "Brier")], row.names = FALSE)

save_table(comparison, "Table4_external_validation.csv")
save_table(data.frame(fold = seq_along(folds), upd_params),
           "TableS_recalibration_parameters.csv")

# ---------------------------------------------------------------------------
# 4. Calibration curves (Figure 4D, 4E)
# ---------------------------------------------------------------------------
calib_curve <- function(p, y, label) {
  # loess-smoothed observed proportion against predicted probability
  df <- data.frame(p = p, y = y)
  fit <- stats::loess(y ~ p, data = df, span = 0.75,
                      control = stats::loess.control(surface = "direct"))
  grid <- seq(min(p), max(p), length.out = 100)
  data.frame(model = label, predicted = grid,
             observed = stats::predict(fit, newdata = data.frame(p = grid)))
}

cal_df <- rbind(calib_curve(pred, y, "Original"),
                calib_curve(pred_upd, y, "Recalibrated"))

p_cal <- ggplot2::ggplot(cal_df, ggplot2::aes(predicted, observed,
                                              colour = model)) +
  ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed",
                       colour = "grey50") +
  ggplot2::geom_line(linewidth = 0.9) +
  ggplot2::coord_cartesian(xlim = c(0, max(c(pred, pred_upd))),
                           ylim = c(0, max(cal_df$observed) * 1.1)) +
  ggplot2::scale_colour_manual(values = c(Original = "grey55",
                                          Recalibrated = "#219EBC")) +
  ggplot2::labs(x = "Predicted probability", y = "Observed probability",
                title = "Calibration: original vs recalibrated") +
  theme_paper()
save_fig(p_cal, "Figure4D_calibration.png", width = 6, height = 5)

# ---------------------------------------------------------------------------
# 5. ROC in the validation set (Figure 4F)
# ---------------------------------------------------------------------------
roc_df <- rbind(
  data.frame(model = sprintf("Original (AUC = %.3f)", as.numeric(pROC::auc(roc_raw))),
             FPR = 1 - roc_raw$specificities, TPR = roc_raw$sensitivities),
  data.frame(model = sprintf("Recalibrated (AUC = %.3f)", as.numeric(pROC::auc(roc_upd))),
             FPR = 1 - roc_upd$specificities, TPR = roc_upd$sensitivities)
)
p_roc <- ggplot2::ggplot(roc_df, ggplot2::aes(FPR, TPR, colour = model)) +
  ggplot2::geom_line(linewidth = 0.9) +
  ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed",
                       colour = "grey50") +
  ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  ggplot2::scale_colour_manual(values = c("grey40", "#219EBC")) +
  ggplot2::labs(x = "1 - Specificity", y = "Sensitivity",
                title = "External validation") +
  theme_paper()
save_fig(p_roc, "Figure4F_roc_external.png", width = 5.5, height = 5)

# ---------------------------------------------------------------------------
# 6. Store predictions for the decision curve
# ---------------------------------------------------------------------------
write_step(list(y = y, original = pred, recalibrated = pred_upd,
                comparison = comparison), "external_predictions")

stamp_session("12_external_validation")
