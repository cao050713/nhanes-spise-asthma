# =============================================================================
# 04_weighted_logistic.R -- survey-weighted association models (Table 2)
#
# Model 1  unadjusted
# Model 2  age, sex, race
# Model 3  Model 2 + fasting glucose, HbA1c, smoking, alcohol, the four
#          leucocyte counts, RDW, platelet count, cardiovascular disease
#          and cancer
#
# Weights are WTMEC2YR / 5 (five pooled cycles); variance estimation uses the
# strata and PSU variables. quasibinomial is used because the survey
# distribution is not a plain binomial likelihood.
#
#   Rscript R/04_weighted_logistic.R
# =============================================================================

source("R/utils.R")
need("survey", "dplyr")

cohort <- read_step("cohort")

options(survey.lonely.psu = "adjust")

design <- survey::svydesign(
  id      = ~SDMVPSU,
  strata  = ~SDMVSTRA,
  weights = ~WTMEC_combind,
  nest    = TRUE,
  data    = cohort
)

f_unadj <- asthma ~ spise
f_demo  <- asthma ~ spise + RIAGENDR + RIDAGEYR + RIDRETH1
f_full  <- asthma ~ spise + RIAGENDR + RIDAGEYR + RIDRETH1 +
  LBXGH + LBXGLU + ALQ110 + SMQ020 +
  LBDMONO + LBDLYMNO + LBDNENO + LBXRDW + LBXPLTSI + LBDEONO +
  CVD + cancer

models <- list(
  `Model 1` = survey::svyglm(f_unadj, family = quasibinomial(), design = design),
  `Model 2` = survey::svyglm(f_demo,  family = quasibinomial(), design = design),
  `Model 3` = survey::svyglm(f_full,  family = quasibinomial(), design = design)
)

# ---------------------------------------------------------------------------
# Odds ratios for SPISE
# ---------------------------------------------------------------------------
spise_or <- do.call(rbind, lapply(names(models), function(m) {
  co <- summary(models[[m]])$coefficients
  b  <- co["spise", "Estimate"]; se <- co["spise", "Std. Error"]
  data.frame(
    Model  = m,
    OR     = exp(b),
    lower  = exp(b - 1.96 * se),
    upper  = exp(b + 1.96 * se),
    p      = co["spise", ncol(co)],
    stringsAsFactors = FALSE
  )
}))
spise_or$OR_95CI <- with(spise_or, fmt_or(OR, lower, upper))
spise_or$p_formatted <- fmt_p(spise_or$p)

print(spise_or[, c("Model", "OR_95CI", "p_formatted")], row.names = FALSE)
save_table(spise_or, "Table2_SPISE_odds_ratios.csv")

# ---------------------------------------------------------------------------
# Full coefficient table for the fully adjusted model
# ---------------------------------------------------------------------------
co <- summary(models$`Model 3`)$coefficients
full_tab <- data.frame(
  term  = rownames(co),
  beta  = co[, "Estimate"],
  se    = co[, "Std. Error"],
  OR    = exp(co[, "Estimate"]),
  lower = exp(co[, "Estimate"] - 1.96 * co[, "Std. Error"]),
  upper = exp(co[, "Estimate"] + 1.96 * co[, "Std. Error"]),
  p     = co[, ncol(co)],
  stringsAsFactors = FALSE
)
full_tab$OR_95CI <- with(full_tab, fmt_or(OR, lower, upper))
full_tab$p_formatted <- fmt_p(full_tab$p)
save_table(full_tab, "Table2_full_model3.csv")

cat("\n--- fully adjusted model ---\n")
print(full_tab[, c("term", "OR_95CI", "p_formatted")], row.names = FALSE)

# ---------------------------------------------------------------------------
# Design-based degrees of freedom, quoted when describing the survey models
# ---------------------------------------------------------------------------
cat("\nDesign df:", survey::degf(design), "\n")

stamp_session("04_weighted_logistic")
