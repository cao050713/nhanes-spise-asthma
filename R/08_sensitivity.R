# =============================================================================
# 08_sensitivity.R -- sensitivity analyses (Supplementary Tables S1 and S2)
#
# Two analyses repeat the same three-model sequence used for Table 2, changing
# one thing at a time:
#
#   (a) SPISE entered as quartiles, Q1 as reference (Table S1)
#   (b) the per-unit model refitted without survey weights (Table S2)
#
# Both are reported for Model 1 (unadjusted), Model 2 (age, sex, race) and
# Model 3 (fully adjusted), because the manuscript describes the gradient
# across all three.
#
#   Rscript R/08_sensitivity.R
# =============================================================================

source("R/utils.R")
need("survey", "dplyr")

cohort <- read_step("cohort")

options(survey.lonely.psu = "adjust")

# ---------------------------------------------------------------------------
# Covariate blocks, identical to R/04
# ---------------------------------------------------------------------------
demo_covars <- c("RIAGENDR", "RIDAGEYR", "RIDRETH1")
full_covars <- c(demo_covars, "LBXGH", "LBXGLU", "ALQ110", "SMQ020",
                 "LBDMONO", "LBDLYMNO", "LBDNENO", "LBXRDW", "LBXPLTSI",
                 "LBDEONO", "CVD", "cancer")

make_design <- function(d) {
  survey::svydesign(id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTMEC_combind,
                    nest = TRUE, data = d)
}

# Pull the rows of a fit whose names start with `prefix` into a tidy frame.
tidy_terms <- function(fit, prefix, model_label, exposure) {
  co <- summary(fit)$coefficients
  rows <- grep(paste0("^", prefix), rownames(co), value = TRUE)
  if (!length(rows)) return(NULL)
  b  <- co[rows, "Estimate"]
  se <- co[rows, "Std. Error"]
  data.frame(
    exposure = exposure,
    Model    = model_label,
    term     = rows,
    OR       = exp(b),
    lower    = exp(b - 1.96 * se),
    upper    = exp(b + 1.96 * se),
    p        = co[rows, ncol(co)],
    stringsAsFactors = FALSE, row.names = NULL
  )
}

# ---------------------------------------------------------------------------
# (a) Quartiles of SPISE, survey-weighted
# ---------------------------------------------------------------------------
# Cut points come from the sample distribution of SPISE. This reproduces the
# published Table S1 (Q4: OR 0.50, p = 0.014). Weighted cut points
# (svyquantile) give a materially different Q4 (0.52, 0.30-0.90); see
# docs/KNOWN_DIFFERENCES.md, B5.
cuts <- unique(stats::quantile(cohort$spise, 0:4 / 4, na.rm = TRUE))
cuts[1] <- cuts[1] - 1e-8
cuts[5] <- cuts[5] + 1e-8

cohort$spise_q <- cut(cohort$spise, breaks = cuts,
                      labels = c("Q1", "Q2", "Q3", "Q4"),
                      include.lowest = TRUE)
cohort$spise_q_num <- as.integer(cohort$spise_q)

cat("Quartile cut points:", paste(round(cuts[2:4], 2), collapse = ", "), "\n")
print(table(cohort$spise_q, useNA = "ifany"))

design_q <- make_design(cohort)

f_q1 <- asthma ~ spise_q
f_q2 <- stats::as.formula(paste("asthma ~ spise_q +", paste(demo_covars, collapse = " + ")))
f_q3 <- stats::as.formula(paste("asthma ~ spise_q +", paste(full_covars, collapse = " + ")))

q_fits <- list(
  `Model 1` = survey::svyglm(f_q1, family = quasibinomial(), design = design_q),
  `Model 2` = survey::svyglm(f_q2, family = quasibinomial(), design = design_q),
  `Model 3` = survey::svyglm(f_q3, family = quasibinomial(), design = design_q)
)

tab_s1 <- do.call(rbind, lapply(names(q_fits),
                                function(m) tidy_terms(q_fits[[m]], "spise_q",
                                                       m, "SPISE quartile")))

# p for trend: quartile rank entered as a single ordinal term
tab_s1_trend <- do.call(rbind, lapply(names(q_fits), function(m) {
  covs <- switch(m,
                 `Model 1` = NULL,
                 `Model 2` = demo_covars,
                 `Model 3` = full_covars)
  f <- stats::as.formula(paste("asthma ~ spise_q_num",
                               if (length(covs)) paste("+", paste(covs, collapse = " + ")) else ""))
  fit <- survey::svyglm(f, family = quasibinomial(), design = design_q)
  co <- summary(fit)$coefficients
  data.frame(exposure = "SPISE quartile", Model = m, term = "p for trend",
             OR = NA_real_, lower = NA_real_, upper = NA_real_,
             p = co["spise_q_num", ncol(co)],
             stringsAsFactors = FALSE, row.names = NULL)
}))

tab_s1 <- rbind(tab_s1, tab_s1_trend)

# ---------------------------------------------------------------------------
# (b) The per-unit model, unweighted
# ---------------------------------------------------------------------------
f_u1 <- asthma ~ spise
f_u2 <- stats::as.formula(paste("asthma ~ spise +", paste(demo_covars, collapse = " + ")))
f_u3 <- stats::as.formula(paste("asthma ~ spise +", paste(full_covars, collapse = " + ")))

u_fits <- list(
  `Model 1` = stats::glm(f_u1, data = cohort, family = binomial()),
  `Model 2` = stats::glm(f_u2, data = cohort, family = binomial()),
  `Model 3` = stats::glm(f_u3, data = cohort, family = binomial())
)

tab_s2 <- do.call(rbind, lapply(names(u_fits),
                                function(m) tidy_terms(u_fits[[m]], "spise",
                                                       m, "SPISE (per unit)")))

# ---------------------------------------------------------------------------
# Formatting and output
# ---------------------------------------------------------------------------
finish <- function(x) {
  x$OR_95CI <- with(x, ifelse(is.na(OR), "-", fmt_or(OR, lower, upper)))
  x$p_formatted <- fmt_p(x$p)
  x[, c("exposure", "Model", "term", "OR_95CI", "p_formatted", "OR", "lower", "upper", "p")]
}

tab_s1 <- finish(tab_s1)
tab_s2 <- finish(tab_s2)

cat("\n--- Table S1. SPISE in quartiles ---\n")
print(tab_s1[, c("Model", "term", "OR_95CI", "p_formatted")], row.names = FALSE)

cat("\n--- Table S2. Per-unit model without survey weights ---\n")
print(tab_s2[, c("Model", "OR_95CI", "p_formatted")], row.names = FALSE)

save_table(tab_s1, "TableS1_sensitivity_quartiles.csv")
save_table(tab_s2, "TableS2_sensitivity_unweighted.csv")

# Size of the shift caused by dropping the weights, for the Results sentence.
w_or <- exp(stats::coef(survey::svyglm(
  stats::as.formula(paste("asthma ~ spise +", paste(full_covars, collapse = " + "))),
  family = quasibinomial(), design = make_design(cohort)))["spise"])
u_or <- exp(stats::coef(u_fits$`Model 3`))["spise"]
cat("\nModel 3, weighted OR:", round(w_or, 3),
    "| unweighted OR:", round(u_or, 3), "\n")

stamp_session("08_sensitivity")
