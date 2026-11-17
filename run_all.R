# =============================================================================
# run_all.R -- run the whole pipeline in order
#
#   Rscript run_all.R                 # download data if missing, then everything
#   Rscript run_all.R --no-download   # skip the download step
#   Rscript run_all.R 02 03 09        # run selected steps only
#
# Intermediate objects are written to data/, so a single step can be re-run
# without repeating the expensive ones (10 and 11 take the longest).
# =============================================================================

STEPS <- c(
  "01_download_data.R",
  "02_build_cohort.R",
  "03_baseline_table.R",
  "04_weighted_logistic.R",
  "05_rcs_breakpoint.R",
  "06_subgroup.R",
  "07_index_comparison.R",
  "08_sensitivity.R",
  "09_feature_selection.R",
  "10_model_comparison_cv.R",
  "11_internal_validation.R",
  "12_external_validation.R",
  "13_decision_curve.R",
  "14_shap.R",
  "15_cv_metrics_recompute.R"
)

args      <- commandArgs(trailingOnly = TRUE)
skip_dl   <- "--no-download" %in% args
requested <- setdiff(args, "--no-download")

steps <- if (skip_dl) setdiff(STEPS, "01_download_data.R") else STEPS

if (length(requested)) {
  ids <- sub("\\.R$", "", sub("^([0-9]+).*$", "\\1", STEPS))
  keep <- ids %in% requested
  steps <- steps[keep]
  if (!length(steps)) {
    stop("No step matched. Available numbers: ", paste(ids, collapse = ", "),
         call. = FALSE)
  }
}

source("R/00_setup.R")

t0 <- Sys.time()
for (s in steps) {
  message("\n================================================================")
  message("Running ", s)
  message("================================================================")
  source(file.path("R", s), echo = FALSE)
}
message(sprintf("\nFinished %d step(s) in %.1f minutes.",
                length(steps),
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))
