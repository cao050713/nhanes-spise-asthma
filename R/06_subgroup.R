# =============================================================================
# 06_subgroup.R -- effect modification (Figure 1C)
#
# Seven interaction tests on the fully adjusted survey model. Age and fasting
# glucose enter the interaction model on their original continuous scale;
# stratum-specific odds ratios are computed for display only.
#
# Interaction p values are judged against 0.05 / 7 = 0.0071.
#
#   Rscript R/06_subgroup.R
# =============================================================================

source("R/utils.R")
need("survey", "dplyr", "ggplot2")

cohort <- read_step("cohort")

options(survey.lonely.psu = "adjust")

design <- survey::svydesign(
  id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTMEC_combind,
  nest = TRUE, data = cohort
)

# Variables that are adjusted for but never used as the stratifier
adj <- c("RIAGENDR", "RIDAGEYR", "RIDRETH1", "SMQ020", "ALQ110", "CVD", "cancer",
         "LBXGLU", "LBXGH", "LBDMONO", "LBDLYMNO", "LBDNENO", "LBXRDW",
        "LBXPLTSI", "LBDEONO")

# ---------------------------------------------------------------------------
# 1. Interaction tests
# ---------------------------------------------------------------------------
interaction_vars <- c("RIAGENDR", "RIDAGEYR", "RIDRETH1", "SMQ020", "ALQ110",
                      "CVD", "cancer")

int_tab <- do.call(rbind, lapply(interaction_vars, function(v) {
  others <- setdiff(adj, v)
  f <- stats::as.formula(sprintf("asthma ~ spise * %s + %s",
                                 v, paste(others, collapse = " + ")))
  fit <- survey::svyglm(f, family = quasibinomial(), design = design)
  p <- survey::regTermTest(fit, stats::as.formula(paste0("~ spise:", v)),
                           method = "Wald")$p
  data.frame(variable = v, p_interaction = p,
             significant = p < ALPHA_SUBGROUP,
             stringsAsFactors = FALSE)
}))
int_tab$p_formatted <- fmt_p(int_tab$p_interaction)

cat("\n--- interaction tests (Bonferroni threshold ",
    round(ALPHA_SUBGROUP, 4), ") ---\n", sep = "")
print(int_tab[, c("variable", "p_formatted", "significant")], row.names = FALSE)
save_table(int_tab, "TableS_interaction_tests.csv")

# ---------------------------------------------------------------------------
# 2. Stratum-specific odds ratios (display only)
# ---------------------------------------------------------------------------
# Each stratifier is cut into display categories. Age and glucose stay
# continuous in the interaction model above; here they are binned only so the
# forest plot has something to show.
cohort <- cohort %>%
  dplyr::mutate(
    age_group = cut(RIDAGEYR, breaks = c(19, 39, 59, Inf),
                    labels = c("20-39", "40-59", "60+")),
    glucose_group = cut(LBXGLU, breaks = c(-Inf, 99.9, 124.9, Inf),
                        labels = c("<100", "100-125", ">=125")),
    race = factor(RIDRETH1, levels = 1:5,
                  labels = c("Mexican American", "Other Hispanic",
                             "Non-Hispanic White", "Non-Hispanic Black",
                             "Other")),
    sex = factor(RIAGENDR, levels = 1:2, labels = c("Male", "Female")),
    smoking = factor(SMQ020, levels = 1:2, labels = c("Never/former", "Current")),
    alcohol = factor(ALQ110, levels = 1:2, labels = c("No", "Yes")),
    cvd = factor(CVD, levels = 0:1, labels = c("No", "Yes")),
    cancer_hist = factor(cancer, levels = 0:1, labels = c("No", "Yes"))
  )

design <- survey::svydesign(
  id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTMEC_combind,
  nest = TRUE, data = cohort
)

strata_map <- list(
  sex     = list(col = "sex",         base = "RIAGENDR"),
  age     = list(col = "age_group",   base = "RIDAGEYR"),
  race    = list(col = "race",        base = "RIDRETH1"),
  smoking = list(col = "smoking",     base = "SMQ020"),
  alcohol = list(col = "alcohol",     base = "ALQ110"),
  cvd     = list(col = "cvd",         base = "CVD"),
  cancer  = list(col = "cancer_hist", base = "cancer")
)

# Within a stratum the model keeps the full adjustment set except the
# stratifying variable itself.
stratum_or <- function(stratum, col, base) {
  keep <- setdiff(adj, base)
  f <- stats::as.formula(paste("asthma ~ spise +", paste(keep, collapse = " + ")))
  lv <- levels(cohort[[col]])

  out <- lapply(lv, function(l) {
    cond  <- call("==", as.name(col), l)
    d_sub <- eval(bquote(subset(design, .(cond))))
    fit <- tryCatch(
      survey::svyglm(f, family = quasibinomial(), design = d_sub),
      error = function(e) NULL
    )
    if (is.null(fit)) return(NULL)
    co <- summary(fit)$coefficients
    if (!"spise" %in% rownames(co)) return(NULL)
    b <- co["spise", "Estimate"]; se <- co["spise", "Std. Error"]
    data.frame(
      stratum = stratum, level = l,
      n = nrow(d_sub$variables), events = sum(d_sub$variables$asthma),
      OR = exp(b), lower = exp(b - 1.96 * se), upper = exp(b + 1.96 * se),
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, Filter(Negate(is.null), out))
}

forest_tab <- do.call(rbind, lapply(names(strata_map), function(k) {
  s <- strata_map[[k]]
  stratum_or(k, s$col, s$base)
}))

print(forest_tab, row.names = FALSE)
save_table(forest_tab, "TableS_subgroup_odds_ratios.csv")

# ---------------------------------------------------------------------------
# 3. Forest plot (Figure 1C)
# ---------------------------------------------------------------------------
forest_tab$label <- paste0(forest_tab$stratum, ": ", forest_tab$level)
forest_tab$label <- factor(forest_tab$label, levels = rev(forest_tab$label))

p <- ggplot2::ggplot(forest_tab,
                     ggplot2::aes(x = OR, y = label, xmin = lower, xmax = upper)) +
  ggplot2::geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
  ggplot2::geom_pointrange(colour = "#E64B35", size = 0.35) +
  ggplot2::labs(x = "Odds ratio per 1-unit higher SPISE (95% CI)", y = NULL,
                title = "Subgroup analyses") +
  theme_paper()

save_fig(p, "Figure1C_subgroups.png", width = 6, height = 6)

stamp_session("06_subgroup")
