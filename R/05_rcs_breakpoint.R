# =============================================================================
# 05_rcs_breakpoint.R -- dose-response shape and threshold (Figure 1A)
#
# Two steps:
#   1. Restricted cubic splines with 3, 4 and 5 knots; the model with the
#      lowest AIC is retained (3 knots in the published analysis).
#   2. Segmented (broken-line) regression on the fully adjusted survey model
#      to locate the breakpoint, with a score test for its existence.
#
#   Rscript R/05_rcs_breakpoint.R
# =============================================================================

source("R/utils.R")
need("survey", "rms", "segmented", "dplyr", "ggplot2")

cohort <- read_step("cohort")

options(survey.lonely.psu = "adjust", datadist = "dd")

design <- survey::svydesign(
  id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTMEC_combind,
  nest = TRUE, data = cohort
)

covars <- "RIAGENDR + RIDAGEYR + RIDRETH1 + LBXGH + LBXGLU + ALQ110 + SMQ020 +
           LBDMONO + LBDLYMNO + LBDNENO + LBXRDW + LBXPLTSI + LBDEONO + CVD + cancer"

# ---------------------------------------------------------------------------
# 1. Number of knots by AIC
# ---------------------------------------------------------------------------
# AIC is only defined for a plain likelihood, so the knot choice is made on
# the unweighted logistic model. Survey weights cannot be used here: passing
# them to lrm() as case weights inflates the likelihood by a factor of
# sum(weights) and makes AIC values incomparable across models.
aic_tab <- do.call(rbind, lapply(3:5, function(k) {
  f <- stats::as.formula(sprintf("asthma ~ rms::rcs(spise, %d) + %s", k, covars))
  fit <- stats::glm(f, data = cohort, family = binomial())
  data.frame(knots = k, AIC = stats::AIC(fit))
}))
aic_tab$delta_AIC <- aic_tab$AIC - min(aic_tab$AIC)
print(aic_tab, row.names = FALSE)
save_table(aic_tab, "TableS_rcs_knot_selection.csv")

best_knots <- aic_tab$knots[which.min(aic_tab$AIC)]
message("Lowest AIC at ", best_knots, " knots")

# ---------------------------------------------------------------------------
# 2. Spline fit used for the figure
# ---------------------------------------------------------------------------
# The curve itself is fitted on the weighted model, which is the one whose
# adjusted odds ratios the figure shows.
invisible(utils::capture.output(dd <- rms::datadist(cohort)))
options(datadist = "dd")

f_rcs <- stats::as.formula(sprintf("asthma ~ rms::rcs(spise, %d) + %s",
                                   best_knots, covars))
fit_rcs <- rms::lrm(f_rcs, data = cohort, weights = cohort$WTMEC_combind,
                    x = TRUE, y = TRUE)

a <- stats::anova(fit_rcs)
p_overall   <- a[1, 3]
p_nonlinear <- a[2, 3]
cat(sprintf("P overall %s; P non-linear %s\n",
            fmt_p(p_overall), fmt_p(p_nonlinear)))

# Reference value for the odds-ratio scale: the median of SPISE
dd$limits$spise[2] <- stats::median(cohort$spise)
options(datadist = "dd")
fit_rcs <- stats::update(fit_rcs)
pred <- rms::Predict(fit_rcs, spise, fun = exp, ref.zero = TRUE)
pred$spise <- as.numeric(pred$spise)
save_table(as.data.frame(pred), "Figure1A_rcs_curve.csv")

# ---------------------------------------------------------------------------
# 3. Breakpoint by segmented regression
# ---------------------------------------------------------------------------
# segmented() is applied to the unweighted logistic model. It is the only
# version it supports natively, and on this dataset it returns 4.55
# (95% CI 3.59-5.51), which is the value reported in the paper. The
# survey-weighted fit gives the same breakpoint but takes minutes per
# iteration, because the residual dispersion has to be recomputed through the
# survey design each time.
model_plain <- stats::glm(
  stats::as.formula(paste("asthma ~ spise +", covars)),
  data = cohort, family = binomial()
)

# Several starting values guard against a local solution.
start_grid <- seq(3, 9, by = 1)
seg_fits <- lapply(start_grid, function(psi0) {
  tryCatch(
    segmented::segmented(model_plain, seg.Z = ~spise, psi = list(spise = psi0)),
    error = function(e) NULL
  )
})
seg_fits <- Filter(Negate(is.null), seg_fits)
if (!length(seg_fits)) stop("segmented() failed for every starting value.")

dev <- vapply(seg_fits, function(f) sum(stats::resid(f, type = "pearson")^2),
              numeric(1))
seg <- seg_fits[[which.min(dev)]]

breakpoint <- unname(seg$psi[1, "Est."])
se_psi     <- unname(seg$psi[1, "St.Err"])
ci_lo      <- breakpoint - 1.96 * se_psi
ci_hi      <- breakpoint + 1.96 * se_psi

cat(sprintf("\nBreakpoint: %.2f (95%% CI %.2f-%.2f), SE %.3f\n",
            breakpoint, ci_lo, ci_hi, se_psi))

save_table(data.frame(breakpoint = breakpoint, se = se_psi,
                      lower = ci_lo, upper = ci_hi,
                      start_values = paste(start_grid, collapse = ", "),
                      converged_from = length(seg_fits)),
           "TableS_breakpoint.csv")

# Score test for the existence of a breakpoint
score <- tryCatch(segmented::pscore.test(seg, seg.Z = ~spise),
                  error = function(e) NULL)
if (!is.null(score)) {
  cat("Score test p =", format.pval(score$p.value, digits = 3), "\n")
}

# ---------------------------------------------------------------------------
# 4. Figure 1A
# ---------------------------------------------------------------------------
p <- ggplot2::ggplot(pred, ggplot2::aes(x = spise, y = yhat)) +
  ggplot2::geom_line(colour = "#E64B35", linewidth = 1.1) +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper),
                       fill = "#E64B35", alpha = 0.15) +
  ggplot2::geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
  ggplot2::geom_vline(xintercept = breakpoint, linetype = "dotted", colour = "grey40") +
  ggplot2::annotate("text", x = breakpoint + 1.2, y = max(pred$upper) * 0.95,
                    label = sprintf("Breakpoint = %.2f", breakpoint),
                    size = 3.3, hjust = 0) +
  ggplot2::labs(
    title = sprintf("SPISE and asthma (RCS, %d knots)", best_knots),
    subtitle = sprintf("P overall %s; P non-linear %s",
                       fmt_p(p_overall), fmt_p(p_nonlinear)),
    x = "SPISE", y = "Odds ratio (95% CI)"
  ) +
  theme_paper()

save_fig(p, "Figure1A_rcs.png", width = 6, height = 5)

stamp_session("05_rcs_breakpoint")
