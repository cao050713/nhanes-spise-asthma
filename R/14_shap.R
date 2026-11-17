# =============================================================================
# 14_shap.R -- SHAP attribution of the locked model (Figure 5)
#
# Kernel SHAP on 100 randomly selected participants from the external
# validation set, with 200 background observations drawn from the
# development set. Attribution is computed on the probability scale from the
# original (pre-recalibration) model; recalibration is a monotone
# transformation of the predicted probability and leaves relative
# contributions unchanged.
#
#   Rscript R/14_shap.R
# =============================================================================

source("R/utils.R")
# caret supplies predict.preProcess for the imputation object stored with the
# locked model; gbm supplies predict.gbm for the model itself.
need("kernelshap", "shapviz", "dplyr", "ggplot2", "gbm", "caret")

split  <- read_step("split")
locked <- read_step("locked_model")

dev <- split$development
val <- split$validation
vars <- locked$vars

set.seed(SEED)

# ---------------------------------------------------------------------------
# 1. Explain a random subset of the validation set
# ---------------------------------------------------------------------------
idx_explain <- sample.int(nrow(val), min(SHAP_N_EXPLAIN, nrow(val)))
idx_bg      <- sample.int(nrow(dev), min(SHAP_N_BACKGROUND, nrow(dev)))

X_explain <- as.data.frame(stats::predict(locked$impute,
                                          val[idx_explain, vars]))
X_bg      <- as.data.frame(stats::predict(locked$impute,
                                          dev[idx_bg, vars]))

# ---------------------------------------------------------------------------
# 2. Kernel SHAP
# ---------------------------------------------------------------------------
pred_fun <- function(object, X) {
  # gbm::predict.gbm wants a data frame whose columns are named like the
  # training data, so X is passed through as-is rather than as a matrix.
  as.numeric(stats::predict(object, as.data.frame(X),
                            n.trees = locked$n_trees, type = "response"))
}

shap <- kernelshap::kernelshap(
  object   = locked$fit,
  X        = X_explain,
  bg_X     = X_bg,
  pred_fun = pred_fun
)

saveRDS(shap, path_data("shap_values.rds"))

# ---------------------------------------------------------------------------
# 3. Relative contribution
# ---------------------------------------------------------------------------
mean_abs <- sort(colMeans(abs(shap$S)), decreasing = TRUE)
share <- data.frame(
  variable   = names(mean_abs),
  mean_abs_SHAP = as.numeric(mean_abs),
  percent    = 100 * as.numeric(mean_abs) / sum(as.numeric(mean_abs))
)
print(share, row.names = FALSE)
save_table(share, "Table5_shap_contribution.csv")

# ---------------------------------------------------------------------------
# 4. Figures
# ---------------------------------------------------------------------------
sv <- shapviz::shapviz(shap, X = X_explain)

p_beeswarm <- shapviz::sv_importance(sv, kind = "beeswarm", max_display = 10) +
  ggplot2::labs(title = "SHAP values (external validation set)") +
  theme_paper()
save_fig(p_beeswarm, "Figure5A_shap_beeswarm.png", width = 6.5, height = 5)

p_bar <- shapviz::sv_importance(sv, kind = "bar", max_display = 10) +
  ggplot2::labs(title = "Mean |SHAP|") +
  theme_paper()
save_fig(p_bar, "Figure5B_shap_bar.png", width = 6, height = 4.5)

if ("spise" %in% colnames(X_explain)) {
  p_dep <- shapviz::sv_dependence(sv, v = "spise", color_var = NULL) +
    ggplot2::labs(title = "SHAP dependence: SPISE") +
    theme_paper()
  save_fig(p_dep, "Figure5C_shap_dependence_spise.png", width = 6, height = 4.5)
}

# ---------------------------------------------------------------------------
# 5. Two individual explanations
# ---------------------------------------------------------------------------
# One participant with a high SPISE and average risk, one with a low SPISE and
# high risk; the paper shows this pair as a waterfall/force plot.
spise_col <- if ("spise" %in% names(X_explain)) "spise" else NULL
if (!is.null(spise_col)) {
  q <- stats::quantile(X_explain[[spise_col]], c(0.1, 0.9), names = FALSE)
  i_high <- which.min(abs(X_explain[[spise_col]] - q[2]))
  i_low  <- which.min(abs(X_explain[[spise_col]] - q[1]))

  for (nm in c(high = i_high, low = i_low)) {
    i <- as.integer(nm)
    p_w <- shapviz::sv_waterfall(sv, row_id = i) +
      ggplot2::labs(title = sprintf("Participant %d (SPISE = %.2f)",
                                    i, X_explain[[spise_col]][i])) +
      theme_paper()
    save_fig(p_w, sprintf("Figure5D_waterfall_%s.png",
                          if (i == i_high) "high_spise" else "low_spise"),
             width = 6.5, height = 4.5)
  }
}

stamp_session("14_shap")
