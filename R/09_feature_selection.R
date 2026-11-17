# =============================================================================
# 09_feature_selection.R -- consensus variable selection (Figure 3)
#
# Three algorithms are applied to the development set only, and only the
# variables confirmed by all three are carried into the models:
#
#   Boruta      all-relevant selection against permuted shadow variables
#   LASSO       glmnet with lambda.min from 10-fold cross-validation
#   RF-RFE      random-forest recursive feature elimination
#
#   Rscript R/09_feature_selection.R
#
# Writes data/features.rds with the three selections and the consensus set.
# =============================================================================

source("R/utils.R")
need("dplyr", "Boruta", "glmnet", "caret", "ranger", "ggplot2", "ggvenn")

split <- read_step("split")
dev   <- split$development

# ---------------------------------------------------------------------------
# 1. Candidate variables
# ---------------------------------------------------------------------------
# Eighteen candidates: the exposure, age, sex, smoking, alcohol intake, the two
# glucose markers, the four leucocyte counts, red cell count, RDW, platelet
# count, total and LDL cholesterol, and the two disease-history items. Race is
# *not* a candidate -- it is an adjustment variable in Models 2 and 3, not a
# predictor the selection step is allowed to keep, and the paper's list of
# twelve retained plus six rejected variables has no race in it.
candidates <- c("spise", "RIAGENDR", "RIDAGEYR", "SMQ020", "ALQ110",
                "CVD", "cancer", "LBXGLU", "LBXGH", "LBDMONO", "LBDLYMNO",
                "LBDNENO", "LBXRBCSI", "LBXRDW", "LBXPLTSI", "LBDEONO",
                "LBXTC", "LBDLDL")

stopifnot(length(candidates) == 18L)

dat <- dev[, c("asthma", candidates)]
dat <- dat[stats::complete.cases(dat), ]
dat$asthma_num <- as.integer(dat$asthma)

cat("Development set:", nrow(dat), "participants,", sum(dat$asthma), "events\n")

# ---------------------------------------------------------------------------
# 2. Boruta
# ---------------------------------------------------------------------------
set.seed(SEED)
bor <- Boruta::Boruta(
  stats::as.formula(paste("asthma ~", paste(candidates, collapse = " + "))),
  data = dat[, c("asthma", candidates)],
  maxRuns = 200, doTrace = 0
)
boruta_feature <- Boruta::getSelectedAttributes(bor, withTentative = FALSE)
cat("\nBoruta:", length(boruta_feature), "variables\n"); print(boruta_feature)

bor_imp <- Boruta::attStats(bor)
bor_imp$variable <- rownames(bor_imp)
bor_imp$decision <- factor(bor_imp$decision,
                           levels = c("Rejected", "Tentative", "Confirmed"))
save_table(bor_imp, "Figure3A_boruta_importance.csv")

p_bor <- ggplot2::ggplot(
  bor_imp[order(bor_imp$medianImp), ],
  ggplot2::aes(x = stats::reorder(variable, medianImp), y = medianImp,
               fill = decision)) +
  ggplot2::geom_col() +
  ggplot2::coord_flip() +
  ggplot2::scale_fill_manual(values = c(Rejected = "grey70",
                                        Tentative = "goldenrod",
                                        Confirmed = "#1B9E77")) +
  ggplot2::labs(x = NULL, y = "Median importance") +
  theme_paper()
save_fig(p_bor, "Figure3A_boruta.png", width = 6, height = 5)

# ---------------------------------------------------------------------------
# 3. LASSO
# ---------------------------------------------------------------------------
set.seed(SEED)
x <- stats::model.matrix(~ . - 1, data = dat[, candidates])
y <- dat$asthma_num

cv_lasso <- glmnet::cv.glmnet(x, y, family = "binomial", alpha = 1,
                              nfolds = 10, standardize = TRUE)
coef_lasso <- stats::coef(cv_lasso, s = "lambda.min")
lasso_feature <- setdiff(rownames(coef_lasso)[as.numeric(coef_lasso) != 0], "(Intercept)")
cat("\nLASSO (lambda.min =", round(cv_lasso$lambda.min, 5), "):",
    length(lasso_feature), "variables\n"); print(lasso_feature)

save_table(data.frame(variable = lasso_feature), "Figure3B_lasso_selected.csv")

p_lasso <- ggplot2::ggplot() +
  ggplot2::geom_point(data = data.frame(lambda = cv_lasso$lambda,
                                        cvm = cv_lasso$cvm),
                      ggplot2::aes(lambda, cvm), size = 0.6, colour = "grey40") +
  ggplot2::geom_vline(xintercept = cv_lasso$lambda.min, linetype = "dashed",
                      colour = "#E64B35") +
  ggplot2::scale_x_log10() +
  ggplot2::labs(x = expression(lambda), y = "Binomial deviance") +
  theme_paper()
save_fig(p_lasso, "Figure3B_lasso_deviance.png", width = 5.5, height = 4.5)

# ---------------------------------------------------------------------------
# 4. Random-forest recursive feature elimination
# ---------------------------------------------------------------------------
set.seed(SEED)
ctrl_rfe <- caret::rfeControl(
  functions = caret::rfFuncs,
  method = "cv", number = 5,
  verbose = FALSE
)
rfe_res <- caret::rfe(
  x = dat[, candidates], y = factor(dat$asthma, levels = c(0, 1),
                                    labels = c("No", "Yes")),
  sizes = 4:16,          # the range searched in the published analysis
  rfeControl = ctrl_rfe
)

# caret always evaluates the full predictor set as well, whatever `sizes`
# says. The published search stopped at 16, so that extra row is dropped
# before the optimum is taken. It matters here: the all-18 model is degenerate
# (it predicts the majority class only, Kappa = 0) and edges out the real
# optimum by 0.0004 of accuracy.
prof <- rfe_res$results
prof_win <- prof[prof$Variables %in% 4:16, ]
best_size <- prof_win$Variables[which.max(prof_win$Accuracy)]
rfe_feature <- rfe_res$variables$var[seq_len(best_size)]

cat("\nRF-RFE: optimum at", best_size, "features\n"); print(rfe_feature)

save_table(prof, "Figure3D_rfe_profile.csv")

p_rfe <- ggplot2::ggplot(prof, ggplot2::aes(Variables, Accuracy)) +
  ggplot2::geom_line(colour = "#E64B35") +
  ggplot2::geom_point(size = 1.6, colour = "#E64B35") +
  ggplot2::geom_vline(xintercept = best_size, linetype = "dashed",
                      colour = "grey50") +
  ggplot2::labs(x = "Number of features", y = "Cross-validated accuracy") +
  theme_paper()
save_fig(p_rfe, "Figure3D_rfe.png", width = 5.5, height = 4.5)

# ---------------------------------------------------------------------------
# 5. Consensus
# ---------------------------------------------------------------------------
reconstructed <- Reduce(intersect, list(boruta_feature, lasso_feature, rfe_feature))
cat("\nReconstructed consensus (", length(reconstructed), " variables ):\n",
    sep = "")
print(reconstructed)

# Boruta reproduces the published set exactly; LASSO and RFE are one variable
# away from it, because both depend on a random fold split whose seed is not
# recorded anywhere. The consensus therefore has to be the archived one --
# everything downstream in R/10-R/14 is built on these six features, and the
# manuscript reports results for these six. See docs/KNOWN_DIFFERENCES.md, A7.
consensus <- CONSENSUS_PUBLISHED

cat("\nConsensus used downstream (", length(consensus), " variables ):\n", sep = "")
print(consensus)
cat("\noverlap with the reconstruction:",
    paste(intersect(consensus, reconstructed), collapse = ", "), "\n")
if (length(setdiff(consensus, reconstructed))) {
  cat("in the archive but not reconstructed:",
      paste(setdiff(consensus, reconstructed), collapse = ", "), "\n")
}
if (length(setdiff(reconstructed, consensus))) {
  cat("reconstructed but not in the archive:",
      paste(setdiff(reconstructed, consensus), collapse = ", "), "\n")
}

sel <- list(boruta = boruta_feature, lasso = lasso_feature,
            rfe = rfe_feature, consensus = consensus,
            consensus_reconstructed = reconstructed)
write_step(sel, "features")

# ---------------------------------------------------------------------------
# 6. Venn diagram (Figure 3E)
# ---------------------------------------------------------------------------
venn_list <- list(Boruta = boruta_feature, LASSO = lasso_feature,
                  `RF-RFE` = rfe_feature)
p_venn <- ggvenn::ggvenn(
  venn_list, show_percentage = FALSE, show_elements = FALSE,
  fill_color = c("#1E90FF", "#4DAF4A", "#E64B35"), fill_alpha = 0.45,
  stroke_color = NA
)
save_fig(p_venn, "Figure3E_venn.png", width = 5, height = 5)

# ---------------------------------------------------------------------------
# 7. Correlation among the selected features (Figure 3F)
# ---------------------------------------------------------------------------
cor_mat <- stats::cor(dev[, consensus], method = "spearman")
print(round(cor_mat, 2))
save_table(as.data.frame(cor_mat), "Figure3F_correlation.csv", row.names = TRUE)

stamp_session("09_feature_selection")
