# =============================================================================
# 11_internal_validation.R -- lock the model, bootstrap optimism (Figure S1)
#
# The algorithm is chosen from the comparison in step 10 by the prespecified
# rule: highest cross-validated AUC, with calibration (Brier, then slope)
# breaking ties. The selected algorithm is then refitted on the whole
# development set *without* resampling, so that its probabilities keep the
# population event rate, and optimism is estimated by bootstrap
# (500 replicates, Harrell's method).
#
#   Rscript R/11_internal_validation.R
# =============================================================================

source("R/utils.R")
need("dplyr", "gbm", "pROC", "caret", "doParallel", "foreach", "ggplot2")

split    <- read_step("split")
features <- read_step("features")
trained  <- read_step("trained_models")

dev  <- split$development
vars <- features$consensus
y    <- as.integer(dev$asthma)

# ---------------------------------------------------------------------------
# 1. Apply the selection rule
# ---------------------------------------------------------------------------
summ <- do.call(rbind, lapply(names(trained), function(nm) {
  p <- trained[[nm]]$pred
  bt <- trained[[nm]]$bestTune
  if (!is.null(bt) && nrow(bt) == 1) {
    keep <- rep(TRUE, nrow(p))
    for (k in names(bt)) {
      if (is.null(p[[k]])) next
      keep <- keep & if (is.numeric(p[[k]])) abs(p[[k]] - bt[[k]]) < 1e-8
                        else as.character(p[[k]]) == as.character(bt[[k]])
    }
    p <- p[keep, , drop = FALSE]
  }
  per <- split(p, p$Resample)
  aucs <- vapply(per, function(d) {
    as.numeric(pROC::auc(pROC::roc(as.integer(d$obs == "Yes"), d$Yes,
                                   quiet = TRUE, direction = "<")))
  }, numeric(1))
  data.frame(Model = nm, AUC = mean(aucs), AUC_sd = stats::sd(aucs))
}))

# Highest AUC wins. Algorithms whose AUC is identical once rounded to three
# decimals count as tied -- that is the grouping the manuscript describes
# ("elastic net, GBM and LASSO were effectively tied" at 0.673, with random
# forest at 0.670 just outside) -- and the Brier score then breaks the tie.
best_auc <- max(summ$AUC)
tied <- summ$Model[round(summ$AUC, 3) == round(best_auc, 3)]
cat("\nBest AUC:", round(best_auc, 3), "| tied algorithms:",
    paste(tied, collapse = ", "), "\n")

if (length(tied) > 1) {
  # Brier on the balanced scale is fine for ranking algorithms trained the
  # same way; see R/15 for the absolute-scale correction.
  brier_by_model <- vapply(tied, function(nm) {
    p <- trained[[nm]]$pred
    bt <- trained[[nm]]$bestTune
    if (!is.null(bt) && nrow(bt) == 1) {
      keep <- rep(TRUE, nrow(p))
      for (k in names(bt)) {
        if (is.null(p[[k]])) next
        keep <- keep & if (is.numeric(p[[k]])) abs(p[[k]] - bt[[k]]) < 1e-8
                          else as.character(p[[k]]) == as.character(bt[[k]])
      }
      p <- p[keep, , drop = FALSE]
    }
    brier(as.integer(p$obs == "Yes"), p$Yes)
  }, numeric(1))
  print(round(brier_by_model, 4))
  selected <- names(which.min(brier_by_model))
} else {
  selected <- tied
}

cat("Selected algorithm:", selected, "\n")
save_table(summ, "Table3_model_selection.csv")

# The bootstrap below uses the gbm interface directly, which is what the
# published model is. If a re-run selects a different algorithm the step has
# to be extended with the matching fit/predict pair; stopping is preferable
# to silently fitting the wrong thing.
if (!identical(selected, "GBM")) {
  stop("The optimism bootstrap in this script is implemented for GBM, the ",
       "algorithm selected in the published analysis, but this run selected '",
       selected, "'.\nExtend fit_final()/predict_final() below for that ",
       "algorithm before continuing.", call. = FALSE)
}

# ---------------------------------------------------------------------------
# 2. Final fit on the development set, without resampling
# ---------------------------------------------------------------------------
best <- trained[[selected]]$bestTune
if (is.null(best)) best <- data.frame()

imp <- caret::preProcess(dev[, vars], method = "medianImpute")
dev_proc <- stats::predict(imp, dev[, vars])

fit_final <- gbm::gbm.fit(
  x = dev_proc, y = y,
  distribution      = "bernoulli",
  n.trees           = if ("n.trees" %in% names(best)) best$n.trees else 200,
  interaction.depth = if ("interaction.depth" %in% names(best)) best$interaction.depth else 2,
  shrinkage         = if ("shrinkage" %in% names(best)) best$shrinkage else 0.05,
  n.minobsinnode    = if ("n.minobsinnode" %in% names(best)) best$n.minobsinnode else 10,
  bag.fraction      = 0.5,
  verbose           = FALSE
)
nt <- if ("n.trees" %in% names(best)) best$n.trees else 200

pred_dev <- stats::predict(fit_final, dev_proc, n.trees = nt, type = "response")

auc_apparent <- as.numeric(pROC::auc(pROC::roc(y, pred_dev, quiet = TRUE,
                                               direction = "<")))
cat(sprintf("\nApparent AUC: %.3f\n", auc_apparent))

# ---------------------------------------------------------------------------
# 3. Bootstrap optimism
# ---------------------------------------------------------------------------
cl <- parallel::makeCluster(max(1, parallel::detectCores() - 1))
doParallel::registerDoParallel(cl)

optimism_vec <- foreach::foreach(
  b = seq_len(N_BOOT_OPTIMISM), .combine = c,
  .packages = c("gbm", "pROC")
) %dopar% {
  set.seed(SEED + b)
  idx <- sample.int(nrow(dev_proc), nrow(dev_proc), replace = TRUE)
  d_b <- dev_proc[idx, , drop = FALSE]
  y_b <- y[idx]
  if (length(unique(y_b)) < 2) return(NA_real_)

  fit_b <- gbm::gbm.fit(
    x = d_b, y = y_b, distribution = "bernoulli",
    n.trees           = nt,
    interaction.depth = if ("interaction.depth" %in% names(best)) best$interaction.depth else 2,
    shrinkage         = if ("shrinkage" %in% names(best)) best$shrinkage else 0.05,
    n.minobsinnode    = if ("n.minobsinnode" %in% names(best)) best$n.minobsinnode else 10,
    bag.fraction = 0.5, verbose = FALSE
  )
  p_boot <- stats::predict(fit_b, d_b, n.trees = nt, type = "response")
  auc_boot <- as.numeric(pROC::auc(pROC::roc(y_b, p_boot, quiet = TRUE,
                                             direction = "<")))
  p_orig <- stats::predict(fit_b, dev_proc, n.trees = nt, type = "response")
  auc_orig <- as.numeric(pROC::auc(pROC::roc(y, p_orig, quiet = TRUE,
                                             direction = "<")))
  auc_boot - auc_orig
}
optimism_vec <- optimism_vec[!is.na(optimism_vec)]

parallel::stopCluster(cl)
foreach::registerDoSEQ()

optimism      <- mean(optimism_vec)
auc_corrected <- auc_apparent - optimism
ci <- stats::quantile(optimism_vec, c(0.025, 0.975))

cat(sprintf("\nBootstrap (%d replicates)\n", length(optimism_vec)))
cat(sprintf("  apparent AUC      : %.3f\n", auc_apparent))
cat(sprintf("  mean optimism     : %.3f (2.5th-97.5th: %.3f to %.3f)\n",
            optimism, ci[1], ci[2]))
cat(sprintf("  optimism-corrected: %.3f\n", auc_corrected))

save_table(data.frame(algorithm = selected, n = nrow(dev),
                      events = sum(y), apparent_AUC = auc_apparent,
                      optimism = optimism, optimism_lower = unname(ci[1]),
                      optimism_upper = unname(ci[2]),
                      corrected_AUC = auc_corrected,
                      n_boot = length(optimism_vec)),
           "Table3_internal_validation.csv")

# ---------------------------------------------------------------------------
# 4. Distribution of the optimism estimates (Figure S1)
# ---------------------------------------------------------------------------
p <- ggplot2::ggplot(data.frame(optimism = optimism_vec),
                     ggplot2::aes(optimism)) +
  ggplot2::geom_histogram(bins = 40, fill = "steelblue", colour = "white") +
  ggplot2::geom_vline(xintercept = optimism, colour = "#E64B35", linewidth = 0.9) +
  ggplot2::geom_vline(xintercept = ci, colour = "#E64B35", linetype = "dashed") +
  ggplot2::labs(x = "Optimism in AUC", y = "Bootstrap replicates",
                title = sprintf("Bootstrap optimism (%d replicates)",
                                length(optimism_vec))) +
  theme_paper()
save_fig(p, "FigureS1_optimism.png", width = 6, height = 4.5)

# ---------------------------------------------------------------------------
# 5. Locked model
# ---------------------------------------------------------------------------
locked <- list(algorithm = selected, fit = fit_final, best = best,
               n_trees = nt, impute = imp, vars = vars)
write_step(locked, "locked_model")

stamp_session("11_internal_validation")
