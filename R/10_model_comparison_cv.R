# =============================================================================
# 10_model_comparison_cv.R -- seven algorithms under repeated cross-validation
#
# Configuration (matches the Methods):
#   repeated 5-fold cross-validation, 10 repetitions, identical fold splits
#   SMOTE applied inside each training fold, held-out fold left untouched
#   seven algorithms: logistic, LASSO, elastic net, random forest,
#                     gradient boosting, k-nearest neighbours, naive Bayes
#
# SMOTE makes the out-of-fold probabilities live on a balanced (50%) prior,
# so the Brier score computed directly from them is not comparable with the
# population event rate. R/15_cv_metrics_recompute.R shifts them back before
# any absolute comparison; within-model ranking is unaffected.
#
#   Rscript R/10_model_comparison_cv.R
# =============================================================================

source("R/utils.R")
need("caret", "dplyr", "pROC", "doParallel", "foreach",
     "ranger", "gbm", "glmnet", "naivebayes", "themis")

split    <- read_step("split")
features <- read_step("features")

dev <- split$development
vars <- features$consensus

dat <- dev[, c("asthma", vars)]
dat$asthma <- factor(dat$asthma, levels = c(0, 1), labels = c("No", "Yes"))

cat("Development set:", nrow(dat), "participants,", sum(dat$asthma == "Yes"),
    "events\n")
cat("Features:", paste(vars, collapse = ", "), "\n")

# ---------------------------------------------------------------------------
# 1. Control object
# ---------------------------------------------------------------------------
# SMOTE is passed as a function rather than as the string "smote". caret's
# built-in string route calls an internal helper (checkInstall) that is not
# exported in current caret builds, so `sampling = "smote"` fails with
# "could not find function checkInstall". The function below is the same
# procedure caret uses internally: themis::smote() on the training fold, with
# the outcome in a temporary column that is removed afterwards.
smote_fold <- function(x, y) {
  dat <- as.data.frame(x)
  dat$.y <- y
  dat <- themis::smote(dat, var = ".y")
  list(x = dat[, !grepl(".y", colnames(dat), fixed = TRUE), drop = FALSE],
       y = dat$.y)
}

# One seed per resample keeps the fold structure identical across algorithms,
# which is what makes the head-to-head comparison fair.
n_resamples <- CV_FOLDS * CV_REPEATS
set.seed(SEED)
seeds <- vector(mode = "list", length = n_resamples + 1)
for (i in seq_len(n_resamples)) seeds[[i]] <- sample.int(1e5, 200)
seeds[[n_resamples + 1]] <- sample.int(1e5, 1)

ctrl <- caret::trainControl(
  method          = "repeatedcv",
  number          = CV_FOLDS,
  repeats         = CV_REPEATS,
  classProbs      = TRUE,
  summaryFunction = caret::twoClassSummary,
  savePredictions = "final",
  sampling        = smote_fold,
  allowParallel   = TRUE,
  seeds           = seeds
)

# ---------------------------------------------------------------------------
# 2. Algorithm grid
# ---------------------------------------------------------------------------
# caret's k-nearest-neighbour wrapper calls knn3() without qualifying it. The
# model specification is loaded from caret's own RData file, so inside a
# parallel worker that call cannot be resolved and caret stops with
# "could not find function knn3". Redefining the fit step with the namespaced
# call is the whole fix.
knn_spec <- caret::getModelInfo("knn")$knn
knn_spec$fit <- function(x, y, wts, param, lev, last, classProbs, ...) {
  caret::knn3(as.matrix(x), y, k = param$k, ...)
}

configs <- list(
  Logistic = list(
    method = "glm", family = "binomial", tuneGrid = NULL
  ),
  LASSO = list(
    method = "glmnet",
    tuneGrid = expand.grid(alpha = 1, lambda = seq(0.001, 0.1, length = 20))
  ),
  ElasticNet = list(
    method = "glmnet",
    tuneGrid = expand.grid(alpha = c(0.2, 0.5, 0.8),
                           lambda = seq(0.001, 0.1, length = 20))
  ),
  RF = list(
    method = "ranger",
    tuneGrid = expand.grid(mtry = c(2, 3, 4, 5),
                           splitrule = "gini",
                           min.node.size = c(1, 5, 10)),
    num.trees = 500, importance = "impurity"
  ),
  GBM = list(
    method = "gbm",
    tuneGrid = expand.grid(n.trees = c(100, 200, 300),
                           interaction.depth = c(1, 2, 3),
                           shrinkage = c(0.01, 0.05, 0.1),
                           n.minobsinnode = c(10, 20)),
    verbose = FALSE
  ),
  KNN = list(
    method = knn_spec,
    tuneGrid = data.frame(k = seq(5, 45, by = 5))
  ),
  NaiveBayes = list(
    method = "naive_bayes",
    tuneGrid = expand.grid(laplace = c(0, 1),
                           usekernel = c(TRUE, FALSE),
                           adjust = c(0.5, 1, 2))
  )
)

# ---------------------------------------------------------------------------
# 3. Fit
# ---------------------------------------------------------------------------
# Parallel workers draw their own random numbers (gbm's bagging subsample,
# for instance), so the default Mersenne-Twister stream would make repeated
# runs differ. L'Ecuyer-CMRG plus a seed before registering the cluster gives
# foreach a reproducible parallel stream; combined with the per-resample seeds
# in the control object the whole step is then deterministic.
RNGkind("L'Ecuyer-CMRG")
set.seed(SEED)

cl <- parallel::makeCluster(max(1, parallel::detectCores() - 1))
doParallel::registerDoParallel(cl)

set.seed(SEED)
trained <- lapply(names(configs), function(nm) {
  message("Training ", nm, " ...")
  cfg <- configs[[nm]]
  args <- list(
    x = dat[, vars], y = dat$asthma,
    method = cfg$method,
    trControl = ctrl,
    metric = "ROC",
    preProcess = c("medianImpute", "center", "scale")
  )
  if (!is.null(cfg$tuneGrid))  args$tuneGrid <- cfg$tuneGrid
  if (!is.null(cfg$family))    args$family <- cfg$family
  if (!is.null(cfg$num.trees)) args$num.trees <- cfg$num.trees
  if (!is.null(cfg$importance)) args$importance <- cfg$importance
  if (!is.null(cfg$verbose))   args$verbose <- cfg$verbose
  fit <- do.call(caret::train, args)
  fit$model_name <- nm
  fit
})
names(trained) <- names(configs)

# Release the workers before the rest of the script runs serially.
parallel::stopCluster(cl)
foreach::registerDoSEQ()

# ---------------------------------------------------------------------------
# 4. Out-of-fold predictions
# ---------------------------------------------------------------------------
oof <- do.call(rbind, lapply(names(trained), function(nm) {
  p <- trained[[nm]]$pred
  bt <- trained[[nm]]$bestTune
  if (!is.null(bt) && nrow(bt) == 1) {
    keep <- rep(TRUE, nrow(p))
    for (k in names(bt)) {
      if (is.null(p[[k]])) next
      keep <- keep & if (is.numeric(p[[k]])) abs(p[[k]] - bt[[k]]) < 1e-8
                        else as.character(p[[k]]) == as.character(bt[[k]])
    }
    p <- p[keep, , drop = FALSE]
  }
  data.frame(Model = nm, Resample = p$Resample, rowIndex = p$rowIndex,
             obs = p$obs, prob = p$Yes, stringsAsFactors = FALSE)
}))

oof$y <- as.integer(oof$obs == "Yes")

# ---------------------------------------------------------------------------
# 5. Headline cross-validated metrics
# ---------------------------------------------------------------------------
metric_one <- function(p, y) {
  p <- clip_prob(p)
  c(AUC = as.numeric(pROC::auc(pROC::roc(y, p, quiet = TRUE, direction = "<"))),
    Brier = brier(y, p),
    Slope = calibration_slope(y, p),
    MeanP = mean(p))
}

per_rep <- oof %>%
  dplyr::group_by(Model, Resample) %>%
  dplyr::group_modify(~ as.data.frame(t(metric_one(.x$prob, .x$y)))) %>%
  dplyr::ungroup()

summ <- per_rep %>%
  dplyr::group_by(Model) %>%
  dplyr::summarise(
    n_repeat = dplyr::n(),
    dplyr::across(c(AUC, Brier, Slope, MeanP),
                  list(mean = mean, sd = stats::sd),
                  .names = "{.col}_{.fn}"),
    .groups = "drop"
  ) %>%
  dplyr::arrange(dplyr::desc(AUC_mean))

cat("\n--- cross-validated performance ---\n")
print(as.data.frame(summ), row.names = FALSE)
save_table(summ, "Table3_cv_performance.csv")
save_table(per_rep, "TableS_cv_by_repeat.csv")

# ---------------------------------------------------------------------------
# 6. Save
# ---------------------------------------------------------------------------
write_step(trained, "trained_models")
write_step(oof, "cv_predictions_oof")

stamp_session("10_model_comparison")
