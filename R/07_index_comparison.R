# =============================================================================
# 07_index_comparison.R -- SPISE vs TyG, HOMA-IR and HOMA-BETA (Figure 1B)
#
# Each index enters the same fully adjusted survey-weighted model, separately.
# Model-based AUCs are compared by bootstrap (500 replicates) with a
# Bonferroni-corrected alpha of 0.05/3 = 0.0167 for the three comparisons
# against SPISE.
#
# The comparison needs fasting insulin. The cohort carries LBXIN, assembled in
# R/02 from GLU_F/GLU_G (2009-2012) and INS_H/I/J (2013-2018); it is missing
# for 55 of 4588 participants (1.2%), so this analysis runs on 98.8% of the
# cohort. The subset size is printed and written to
# output/tables/index_comparison_sample.csv.
#
#   Rscript R/07_index_comparison.R
# =============================================================================

source("R/utils.R")
need("dplyr", "survey", "pROC", "ggplot2")

cohort <- read_step("cohort")

# ---------------------------------------------------------------------------
# 1. HOMA indices
# ---------------------------------------------------------------------------
dat <- cohort %>%
  dplyr::filter(!is.na(LBXIN), LBXGLU > 63) %>%      # HOMA-BETA denominator
  dplyr::mutate(
    HOMA_IR   = LBXIN * LBXGLU / 405,
    HOMA_BETA = 360 * LBXIN / (LBXGLU - 63),
    TyG       = log(LBXTR * LBXGLU / 2)
  )

sample_tab <- data.frame(
  n_full_cohort = nrow(cohort),
  n_with_insulin = nrow(dat),
  events = sum(dat$asthma),
  cycles = paste(sort(unique(dat$cycle)), collapse = "; "),
  stringsAsFactors = FALSE
)
print(sample_tab, row.names = FALSE)
save_table(sample_tab, "index_comparison_sample.csv")

# ---------------------------------------------------------------------------
# 2. One survey model per index
# ---------------------------------------------------------------------------
options(survey.lonely.psu = "adjust")
design <- survey::svydesign(
  id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTMEC_combind,
  nest = TRUE, data = dat
)

covars <- "RIAGENDR + RIDAGEYR + RIDRETH1 + LBXGH + LBXGLU + ALQ110 + SMQ020 +
           LBDMONO + LBDLYMNO + LBDNENO + LBXRDW + LBXPLTSI + LBDEONO + CVD + cancer"

indices <- c("spise", "TyG", "HOMA_IR", "HOMA_BETA")
labels  <- c(spise = "SPISE", TyG = "TyG", HOMA_IR = "HOMA-IR",
             HOMA_BETA = "HOMA-BETA")

preds <- lapply(indices, function(v) {
  f <- stats::as.formula(sprintf("asthma ~ %s + %s", v, covars))
  fit <- survey::svyglm(f, family = quasibinomial(), design = design)
  as.numeric(stats::predict(fit, type = "response"))
})
names(preds) <- labels
y <- dat$asthma

auc_tab <- do.call(rbind, lapply(names(preds), function(k) {
  r <- pROC::roc(y, preds[[k]], quiet = TRUE, direction = "<")
  ci <- as.numeric(pROC::ci.auc(r))
  data.frame(
    index = k,
    AUC = as.numeric(pROC::auc(r)),
    lower = ci[1], upper = ci[3],
    stringsAsFactors = FALSE
  )
}))
auc_tab$AUC_95CI <- with(auc_tab,
                         sprintf("%.3f (%.3f-%.3f)", AUC, lower, upper))
print(auc_tab, row.names = FALSE)
save_table(auc_tab, "Table2_index_AUCs.csv")

# ---------------------------------------------------------------------------
# 3. Bootstrap comparison against SPISE
# ---------------------------------------------------------------------------
set.seed(SEED)

boot_delta <- function(p_ref, p_alt, y, B = N_BOOT_AUC) {
  n <- length(y)
  obs <- as.numeric(pROC::auc(pROC::roc(y, p_ref, quiet = TRUE, direction = "<"))) -
         as.numeric(pROC::auc(pROC::roc(y, p_alt, quiet = TRUE, direction = "<")))
  d <- numeric(B)
  for (b in seq_len(B)) {
    i <- sample.int(n, n, replace = TRUE)
    if (length(unique(y[i])) < 2) { d[b] <- NA; next }
    a1 <- as.numeric(pROC::auc(pROC::roc(y[i], p_ref[i], quiet = TRUE, direction = "<")))
    a2 <- as.numeric(pROC::auc(pROC::roc(y[i], p_alt[i], quiet = TRUE, direction = "<")))
    d[b] <- a1 - a2
  }
  d <- d[!is.na(d)]
  # two-sided bootstrap p: proportion of replicates on the other side of zero
  p <- 2 * min(mean(d <= 0), mean(d >= 0))
  list(delta = obs,
       lower = unname(stats::quantile(d, 0.025)),
       upper = unname(stats::quantile(d, 0.975)),
       p = min(p, 1))
}

cmp <- do.call(rbind, lapply(labels[-1], function(k) {
  r <- boot_delta(preds[["SPISE"]], preds[[k]], y)
  data.frame(
    comparison = paste0("SPISE vs ", k),
    delta_AUC = r$delta, lower = r$lower, upper = r$upper, p = r$p,
    stringsAsFactors = FALSE
  )
}))
cmp$p_bonferroni <- pmin(cmp$p * 3, 1)
cmp$significant <- cmp$p_bonferroni < 0.05
cmp$delta_95CI <- with(cmp, sprintf("%.3f (%.3f to %.3f)", delta_AUC, lower, upper))
cmp$p_formatted <- fmt_p(cmp$p_bonferroni)

cat("\n--- AUC differences vs SPISE (", N_BOOT_AUC, " bootstrap replicates, ",
    "Bonferroni-adjusted) ---\n", sep = "")
print(cmp[, c("comparison", "delta_95CI", "p_formatted", "significant")],
      row.names = FALSE)
save_table(cmp, "Table2_index_AUC_comparisons.csv")

# ---------------------------------------------------------------------------
# 4. ROC curves (Figure 1B)
# ---------------------------------------------------------------------------
roc_df <- do.call(rbind, lapply(names(preds), function(k) {
  r <- pROC::roc(y, preds[[k]], quiet = TRUE, direction = "<")
  a <- auc_tab$AUC[auc_tab$index == k]
  data.frame(
    index = sprintf("%s (AUC = %.3f)", k, a),
    FPR = 1 - r$specificities,
    TPR = r$sensitivities,
    stringsAsFactors = FALSE
  )
}))
order_lv <- sprintf("%s (AUC = %.3f)", auc_tab$index, auc_tab$AUC)
roc_df$index <- factor(roc_df$index, levels = order_lv)

p <- ggplot2::ggplot(roc_df, ggplot2::aes(FPR, TPR, colour = index)) +
  ggplot2::geom_line(linewidth = 0.8) +
  ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed",
                       colour = "grey50") +
  ggplot2::coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  ggplot2::scale_colour_manual(
    values = stats::setNames(unname(PAL_MODEL[c("SPISE", "TyG", "HOMA_IR",
                                                "HOMA_BETA")]), order_lv)
  ) +
  ggplot2::labs(x = "1 - Specificity", y = "Sensitivity",
                title = "Discrimination of four insulin-sensitivity indices") +
  theme_paper() +
  ggplot2::theme(legend.position = c(0.78, 0.18),
                 legend.background = ggplot2::element_rect(
                   fill = "white", colour = "grey80"))

save_fig(p, "Figure1B_roc_indices.png", width = 6, height = 5.5)

stamp_session("07_index_comparison")
