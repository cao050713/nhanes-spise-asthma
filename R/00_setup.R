# =============================================================================
# 00_setup.R -- environment check
#
# Run this first (or just run `source("R/00_setup.R")` in an interactive
# session). It verifies that every package used by the pipeline is installed
# and prints the settings that the analysis will use.
#
#   Rscript R/00_setup.R
# =============================================================================

source("R/utils.R")

pkgs <- c(
  "haven",         # read NHANES .XPT files
  "dplyr",         # data manipulation
  "survey",        # survey-weighted models
  "segmented",     # breakpoint estimation
  "pROC",          # ROC / AUC
  "caret",         # cross-validation framework
  "glmnet",        # LASSO and elastic net
  "ranger",        # random forest (caret method "ranger")
  "randomForest",  # random forest (Boruta, RF-RFE)
  "gbm",           # gradient boosting
  "naivebayes",    # naive Bayes (caret method "naive_bayes")
  "Boruta",        # all-relevant feature selection
  "kernelshap",    # SHAP values
  "shapviz",       # SHAP figures
  "ggplot2",       # figures
  "ggvenn",        # Venn diagram
  "doParallel", "foreach"
)

need(pkgs)

cat("\n--- environment ---\n")
cat("R version      :", R.version.string, "\n")
cat("project root   :", PROJ_ROOT, "\n")

cat("\n--- analysis settings (must match the Methods) ---\n")
cat("seed                     :", SEED, "\n")
cat("cycles                   :", N_CYCLES, "(weight divisor", WEIGHT_DIVISOR, ")\n")
cat("age >=                   :", AGE_MIN, "\n")
cat("development set          : cycles ending", DEV_LAST_YEAR, "\n")
cat("validation set           :", paste(VAL_YEARS, collapse = "-"), "\n")
cat("CV                       :", CV_FOLDS, "folds x", CV_REPEATS, "repetitions\n")
cat("bootstrap (optimism)     :", N_BOOT_OPTIMISM, "\n")
cat("bootstrap (AUC compare)  :", N_BOOT_AUC, "\n")
cat("Bonferroni alpha (index) :", round(ALPHA_INDEX, 4), "\n")
cat("Bonferroni alpha (subgr) :", round(ALPHA_SUBGROUP, 4), "\n")
cat("SHAP explained / background:", SHAP_N_EXPLAIN, "/", SHAP_N_BACKGROUND, "\n")
cat("RCS knots                :", RCS_KNOTS, "\n")

cat("\n--- data ---\n")
raw_files <- list.files(path_raw(), pattern = "\\.XPT$", ignore.case = TRUE)
cat("raw XPT files in data/raw:", length(raw_files), "\n")
if (!length(raw_files)) {
  cat("None yet -- run R/01_download_data.R to fetch them from the NHANES website.\n")
}
cohort_file <- path_data("cohort.rds")
cat("derived cohort:", if (file.exists(cohort_file)) "present" else "not built yet",
    "\n")

stamp_session("00_setup")
cat("\nSetup OK.\n")
