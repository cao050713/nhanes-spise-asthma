# =============================================================================
# utils.R -- shared helpers and global constants
#
# Sourced by every analysis script. Defines the project root, the analysis
# constants that must match the published Methods, and small helpers for
# reading/writing intermediate objects and formatting results.
#
# Nothing in this file runs an analysis on its own.
# =============================================================================

# ---------------------------------------------------------------------------
# 1. Project root
# ---------------------------------------------------------------------------
# Scripts are meant to be run from the repository root (see README).
# If that is not the case we walk up from the current directory looking for
# the marker file, so `Rscript R/03_baseline_table.R` also works from R/.
.find_root <- function(start = getwd()) {
  d <- normalizePath(start, mustWork = FALSE)
  for (i in seq_len(6)) {
    if (file.exists(file.path(d, ".project-root"))) return(d)
    parent <- dirname(d)
    if (identical(parent, d)) break
    d <- parent
  }
  stop("Could not locate project root (no .project-root marker found). ",
       "Run the script from the repository root, or run R/00_setup.R first.")
}

PROJ_ROOT <- .find_root()

path_data   <- function(...) file.path(PROJ_ROOT, "data", ...)
path_out    <- function(...) file.path(PROJ_ROOT, "output", ...)
path_fig    <- function(...) file.path(PROJ_ROOT, "output", "figures", ...)
path_tab    <- function(...) file.path(PROJ_ROOT, "output", "tables", ...)
path_raw    <- function(...) file.path(PROJ_ROOT, "data", "raw", ...)

for (.d in c(path_data(), path_fig(), path_tab(), path_raw())) {
  dir.create(.d, showWarnings = FALSE, recursive = TRUE)
}
rm(.d)

# ---------------------------------------------------------------------------
# 2. Analysis constants
# ---------------------------------------------------------------------------
# Every number here is one that also appears in the Methods section of the
# paper. Change them only together with the manuscript.

SEED            <- 42L     # master seed for every stochastic step

# Survey design
N_CYCLES        <- 5L      # 2009-2010 ... 2017-2018
WEIGHT_DIVISOR  <- N_CYCLES

# Cohort definition
AGE_MIN         <- 20L
DEV_LAST_YEAR   <- 2016L   # development set = cycles ending 2009-2016
VAL_YEARS       <- c(2017L, 2018L)

# Cross-validation of the algorithm comparison
CV_FOLDS        <- 5L
CV_REPEATS      <- 10L

# Bootstrap
N_BOOT_OPTIMISM <- 500L    # internal validation (optimism correction)
N_BOOT_AUC      <- 500L    # AUC comparison across the four indices

# Multiplicity
ALPHA_INDEX     <- 0.05 / 3    # SPISE vs TyG, HOMA-IR, HOMA-BETA
ALPHA_SUBGROUP  <- 0.05 / 7    # seven interaction tests

# Decision curve analysis
DCA_THRESHOLDS  <- seq(0.02, 0.40, by = 0.01)
DCA_REPORTED    <- c(0.10, 0.15)

# SHAP
SHAP_N_EXPLAIN  <- 100L    # participants explained in the external set
SHAP_N_BACKGROUND <- 200L

# RCS
RCS_KNOTS       <- 3L

# The six features the published model was built on. The three selection
# algorithms are stochastic and the original run's seeds were never recorded,
# so R/09 reconstructs them as a check but the archived set below is what the
# models in R/10-R/14 use. See docs/KNOWN_DIFFERENCES.md, item A7.
CONSENSUS_PUBLISHED <- c("spise", "LBDLYMNO", "LBDNENO", "LBDEONO",
                         "cancer", "CVD")

# Same reasoning for the algorithm that was locked: R/11 applies the
# prespecified selection rule and reports its winner, but fits and locks the
# algorithm the manuscript reports. See docs/KNOWN_DIFFERENCES.md, item A7.
ALGORITHM_PUBLISHED <- "GBM"

# ---------------------------------------------------------------------------
# 3. Intermediate objects
# ---------------------------------------------------------------------------
# Each step reads the object written by the previous one, so a step can be
# re-run without redoing everything upstream.

read_step <- function(name) {
  f <- path_data(paste0(name, ".rds"))
  if (!file.exists(f)) {
    stop("Missing intermediate object: ", f, "\n",
         "Run the upstream script first (see README for the order).",
         call. = FALSE)
  }
  readRDS(f)
}

write_step <- function(object, name) {
  saveRDS(object, path_data(paste0(name, ".rds")))
  invisible(object)
}

# ---------------------------------------------------------------------------
# 4. Packages
# ---------------------------------------------------------------------------
# Called at the top of every script: fails early with a readable message
# instead of a cryptic "could not find function".
need <- function(...) {
  pkgs <- c(...)
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) {
    stop("Missing package(s): ", paste(missing, collapse = ", "),
         "\nInstall with: install.packages(c(",
         paste0('"', missing, '"', collapse = ", "), "))",
         call. = FALSE)
  }
  invisible(TRUE)
}

`%||%` <- function(x, y) if (is.null(x)) y else x

# dplyr's pipe, re-exported so the scripts can use it without attaching the
# package (every dplyr call is still written with the dplyr:: prefix).
if (requireNamespace("dplyr", quietly = TRUE)) {
  `%>%` <- dplyr::`%>%`
}

# Same reason: `%dopar%` lives in foreach, and the parallel scripts call
# foreach::foreach() without attaching the package.
if (requireNamespace("foreach", quietly = TRUE)) {
  `%dopar%` <- foreach::`%dopar%`
}

# ---------------------------------------------------------------------------
# 5. Formatting helpers
# ---------------------------------------------------------------------------
fmt_p <- function(p, digits = 3) {
  ifelse(is.na(p), NA_character_,
         ifelse(p < 0.001, "< 0.001",
                paste0("= ", formatC(p, format = "f", digits = digits))))
}

fmt_or <- function(est, lo, hi, digits = 2) {
  paste0(formatC(est, format = "f", digits = digits), " (",
         formatC(lo,  format = "f", digits = digits), "-",
         formatC(hi,  format = "f", digits = digits), ")")
}

# ---------------------------------------------------------------------------
# 6. Prediction-metric helpers
# ---------------------------------------------------------------------------
clip_prob <- function(p, eps = 1e-6) pmin(pmax(p, eps), 1 - eps)
logit     <- function(p) qlogis(clip_prob(p))

# Calibration-in-the-large and slope of a set of predicted probabilities.
calibration_slope <- function(y, p) {
  fit <- suppressWarnings(glm(y ~ logit(p), family = binomial()))
  unname(coef(fit)[2])
}

calibration_intercept <- function(y, p) {
  fit <- suppressWarnings(glm(y ~ logit(p), family = binomial()))
  unname(coef(fit)[1])
}

brier <- function(y, p) mean((p - y)^2)

# Brier score of the intercept-only model, i.e. the value a model with no
# information would obtain at this event rate.
brier_null <- function(y) mean(y) * (1 - mean(y))

# SMOTE prior shift.
#
# When SMOTE is applied inside the resampling loop the training folds are
# balanced to 1:1, so the out-of-fold probabilities live on a 50% prior. To
# compare them with a Brier score computed at the population event rate, the
# log-odds have to be moved back by logit(prev_population) - logit(0.5).
prior_shift <- function(prev, prev_train = 0.5) qlogis(prev) - qlogis(prev_train)

# ---------------------------------------------------------------------------
# 7. Plot defaults
# ---------------------------------------------------------------------------
PAL_MODEL <- c(
  SPISE     = "#0072B2",
  TyG       = "#E69F00",
  HOMA_IR   = "#009E73",
  HOMA_BETA = "#D55E00"
)

theme_paper <- function(base_size = 11) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      plot.title   = ggplot2::element_text(face = "bold", hjust = 0.5),
      axis.title   = ggplot2::element_text(face = "bold"),
      legend.title = ggplot2::element_blank()
    )
}

save_fig <- function(plot, filename, width = 6, height = 5) {
  ggplot2::ggsave(path_fig(filename), plot, width = width, height = height,
                  dpi = 300, bg = "white")
  ggplot2::ggsave(sub("\\.png$", ".pdf", path_fig(filename)), plot,
                  width = width, height = height, bg = "white")
  invisible(filename)
}

save_table <- function(x, filename, row.names = FALSE) {
  utils::write.csv(x, path_tab(filename), row.names = row.names)
  invisible(filename)
}

# ---------------------------------------------------------------------------
# 8. Session stamp
# ---------------------------------------------------------------------------
stamp_session <- function(label) {
  f <- file.path(PROJ_ROOT, "logs", paste0(label, "_sessionInfo.txt"))
  utils::capture.output(utils::sessionInfo(), file = f)
  invisible(f)
}
