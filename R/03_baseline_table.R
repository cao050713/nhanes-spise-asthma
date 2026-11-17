# =============================================================================
# 03_baseline_table.R -- Table 1
#
# Baseline characteristics by asthma status in the full analytic sample.
# Continuous variables are summarised as median (IQR) and compared with the
# Mann-Whitney U test; categorical variables as n (%) and compared with the
# chi-square test, or Fisher's exact test when expected counts are small.
# No survey weights are used here, matching the published Table 1.
#
#   Rscript R/03_baseline_table.R
# =============================================================================

source("R/utils.R")
need("dplyr")

cohort <- read_step("cohort")

cont_vars <- c(
  age            = "RIDAGEYR",
  bmi            = "BMXBMI",
  spise          = "spise",
  fasting_glucose = "LBXGLU",
  hba1c          = "LBXGH",
  monocyte       = "LBDMONO",
  lymphocyte     = "LBDLYMNO",
  neutrophil     = "LBDNENO",
  rdw            = "LBXRDW",
  platelet       = "LBXPLTSI",
  eosinophil     = "LBDEONO",
  total_cholesterol = "LBXTC",
  triglycerides  = "LBXTR",
  hdl_cholesterol = "LBDHDD",
  ldl_cholesterol = "LBDLDL"
)

cat_vars <- c(
  sex      = "RIAGENDR",
  race     = "RIDRETH1",
  smoking  = "SMQ020",
  alcohol  = "ALQ110",
  cvd      = "CVD",
  cancer   = "cancer"
)

describe_continuous <- function(var, label) {
  x <- cohort[[var]]; y <- cohort$asthma
  g0 <- x[y == 0]; g1 <- x[y == 1]

  q <- function(v) {
    if (!length(v)) return(c(NA, NA))
    stats::quantile(v, c(0.25, 0.5, 0.75), na.rm = TRUE, names = FALSE)
  }
  q0 <- q(g0); q1 <- q(g1); qa <- q(x)   # qa[2] is the median

  p <- suppressWarnings(stats::wilcox.test(g0, g1, exact = FALSE)$p.value)

  data.frame(
    Characteristic = label,
    Level = "median (IQR)",
    No_asthma = sprintf("%.2f (%.2f-%.2f)", q0[2], q0[1], q0[3]),
    Asthma    = sprintf("%.2f (%.2f-%.2f)", q1[2], q1[1], q1[3]),
    Total     = sprintf("%.2f (%.2f-%.2f)", qa[2], qa[1], qa[3]),
    p_value   = p,
    stringsAsFactors = FALSE
  )
}

describe_categorical <- function(var, label) {
  x <- cohort[[var]]; y <- cohort$asthma
  lv <- sort(unique(stats::na.omit(x)))

  rows <- lapply(lv, function(l) {
    n0 <- sum(x == l & y == 0); n1 <- sum(x == l & y == 1)
    d0 <- sum(y == 0);          d1 <- sum(y == 1)
    data.frame(
      Characteristic = label,
      Level = as.character(l),
      No_asthma = sprintf("%d (%.1f%%)", n0, 100 * n0 / d0),
      Asthma    = sprintf("%d (%.1f%%)", n1, 100 * n1 / d1),
      Total     = sprintf("%d (%.1f%%)", n0 + n1,
                          100 * (n0 + n1) / length(y)),
      p_value   = NA_real_,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)

  tbl <- table(x, y)
  expected <- suppressWarnings(stats::chisq.test(tbl)$expected)
  use_fisher <- any(expected < 1) ||
    mean(expected < 5) > 0.2
  p <- if (use_fisher) {
    stats::fisher.test(tbl)$p.value
  } else {
    suppressWarnings(stats::chisq.test(tbl)$p.value)
  }
  out$p_value[1] <- p
  out$test <- if (use_fisher) "Fisher exact" else "Chi-square"
  out
}

cont_tab <- do.call(rbind, lapply(names(cont_vars),
                                  function(k) describe_continuous(cont_vars[[k]], k)))
cont_tab$test <- "Mann-Whitney U"

cat_tab <- do.call(rbind, lapply(names(cat_vars),
                                 function(k) describe_categorical(cat_vars[[k]], k)))

table1 <- rbind(cont_tab, cat_tab)
table1$p_value_formatted <- fmt_p(table1$p_value)

print(table1, row.names = FALSE)
save_table(table1, "Table1_baseline.csv")

# The Results prose quotes means, not medians, when it describes age, BMI and
# SPISE ("33.2 vs 29.4 kg/m2", "5.20 vs 6.03"); the published table itself
# shows medians. Both summaries are written out so either can be checked.
# See docs/KNOWN_DIFFERENCES.md, B7.
mean_tab <- do.call(rbind, lapply(names(cont_vars), function(k) {
  x <- cohort[[cont_vars[[k]]]]; y <- cohort$asthma
  m <- function(v) sprintf("%.2f (%.2f)", mean(v, na.rm = TRUE), stats::sd(v, na.rm = TRUE))
  p <- suppressWarnings(stats::wilcox.test(x[y == 0], x[y == 1], exact = FALSE)$p.value)
  data.frame(Characteristic = k, Level = "mean (SD)",
             No_asthma = m(x[y == 0]), Asthma = m(x[y == 1]),
             Total = m(x), p_value = p, test = "Mann-Whitney U",
             stringsAsFactors = FALSE)
}))
mean_tab$p_value_formatted <- fmt_p(mean_tab$p_value)
print(mean_tab[, c("Characteristic", "No_asthma", "Asthma", "p_value_formatted")],
      row.names = FALSE)
save_table(mean_tab, "Table1_baseline_means.csv")

# Overall counts, as reported in the first line of the Results
cat(sprintf("\nTotal n = %d; asthma n = %d (%.1f%%)\n",
            nrow(cohort), sum(cohort$asthma),
            100 * mean(cohort$asthma)))

stamp_session("03_baseline")
