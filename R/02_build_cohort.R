# =============================================================================
# 02_build_cohort.R -- assemble the analytic cohort
#
# Reads the raw NHANES .XPT files, derives SPISE, the asthma outcome and the
# covariates, applies the eligibility rules, and writes
#
#     data/cohort.rds          one row per participant (all five cycles)
#     data/split.rds           the development / validation split
#
#   Rscript R/02_build_cohort.R
#
# Expected result with the public 2009-2018 releases:
#     n = 4588 participants, 413 with asthma
#     development set (2009-2016) n = 2445, 205 events
#     validation set  (2017-2018) n = 2143, 208 events
# The script checks these numbers and stops if they do not match, because a
# mismatch means the input files or the eligibility rules differ from the
# ones used in the paper.
# =============================================================================

source("R/utils.R")
need("haven", "dplyr")

# ---------------------------------------------------------------------------
# 1. Asthma definition
# ---------------------------------------------------------------------------
# MCQ010  "Has a doctor or other health professional ever told you that you
#          had asthma?"
# MCQ035  "Do you still have asthma?"   (asked only when MCQ010 = 1)
#
# The published analysis uses *current* asthma: a positive answer to both
# items. This gives 413 cases in the analytic sample. Using MCQ010 alone
# gives 637 cases and does not reproduce the paper. See
# docs/KNOWN_DIFFERENCES.md, item 1.
asthma_from_items <- function(mcq010, mcq035) {
  dplyr::case_when(
    mcq010 == 1 & mcq035 == 1 ~ 1L,   # ever diagnosed and still has asthma
    mcq010 == 1 & mcq035 == 2 ~ 0L,   # ever diagnosed, no longer has it
    mcq010 == 2               ~ 0L,   # never diagnosed
    TRUE                      ~ NA_integer_
  )
}

cycles <- c(F = "2009-2010", G = "2011-2012", H = "2013-2014",
            I = "2015-2016", J = "2017-2018")

# ---------------------------------------------------------------------------
# 2. Read one cycle
# ---------------------------------------------------------------------------
read_cycle <- function(suffix) {

  rd <- function(component, vars) {
    f <- path_raw(sprintf("%s_%s.XPT", component, suffix))
    if (!file.exists(f)) {
      stop("Missing raw file: ", f, "\nRun R/01_download_data.R first.",
           call. = FALSE)
    }
    x <- haven::read_xpt(f)
    keep <- intersect(vars, names(x))
    if (!length(keep)) {
      stop("None of the expected variables found in ", basename(f),
           "\nFound: ", paste(names(x), collapse = ", "), call. = FALSE)
    }
    x[, keep, drop = FALSE]
  }

  demo <- rd("DEMO", c("SEQN", "RIAGENDR", "RIDAGEYR", "RIDRETH1",
                       "SDMVPSU", "SDMVSTRA", "WTMEC2YR"))
  bmx  <- rd("BMX",  c("SEQN", "BMXBMI"))
  hdl  <- rd("HDL",  c("SEQN", "LBDHDD"))
  tri  <- rd("TRIGLY", c("SEQN", "LBXTR", "LBDLDL"))
  tch  <- rd("TCHOL", c("SEQN", "LBXTC"))
  glu_raw <- haven::read_xpt(path_raw(sprintf("GLU_%s.XPT", suffix)))
  glu_keep <- intersect(c("SEQN", "LBXGLU", "LBXIN"), names(glu_raw))
  if (!"LBXGLU" %in% glu_keep) {
    stop("LBXGLU not found in ", basename(sprintf("GLU_%s.XPT", suffix)),
         "\nFound: ", paste(names(glu_raw), collapse = ", "), call. = FALSE)
  }
  glu <- glu_raw[, glu_keep, drop = FALSE]

  # Fasting insulin: part of the glucose file up to 2011-2012, its own INS
  # component from 2013-2014 onward (which is why R/01 downloads INS_H/I/J).
  if ("LBXIN" %in% names(glu)) {
    ins <- glu[, c("SEQN", "LBXIN"), drop = FALSE]
    glu <- glu[, c("SEQN", "LBXGLU"), drop = FALSE]
  } else {
    ins <- rd("INS", c("SEQN", "LBXIN"))
  }

  ghb  <- rd("GHB",  c("SEQN", "LBXGH"))
  cbc  <- rd("CBC",  c("SEQN", "LBDMONO", "LBDLYMNO", "LBDNENO", "LBXRBCSI",
                       "LBXRDW", "LBXPLTSI", "LBDEONO"))
  smq  <- rd("SMQ",  c("SEQN", "SMQ020"))
  mcq  <- rd("MCQ",  c("SEQN", "MCQ010", "MCQ035", "MCQ160B", "MCQ160C",
                       "MCQ160D", "MCQ160F", "MCQ220"))

  # Alcohol item was renamed between cycles: ALQ110 (2009-2016) -> ALQ111 (2017-2018)
  alq_f   <- path_raw(sprintf("ALQ_%s.XPT", suffix))
  alq_raw <- haven::read_xpt(alq_f)
  alq_var <- if ("ALQ110" %in% names(alq_raw)) "ALQ110" else "ALQ111"
  alq     <- alq_raw[, c("SEQN", alq_var), drop = FALSE]
  if (alq_var == "ALQ111") names(alq)[2] <- "ALQ110"

  out <- Reduce(function(a, b) dplyr::inner_join(a, b, by = "SEQN"),
                list(demo, bmx, hdl, tri, tch, glu, ins, ghb, cbc, smq, alq,
                     mcq))
  out$cycle <- cycles[[suffix]]
  out
}

# ---------------------------------------------------------------------------
# 3. Derive variables and apply eligibility
# ---------------------------------------------------------------------------
raw <- dplyr::bind_rows(lapply(names(cycles), read_cycle))
message("Rows after merging all components: ", nrow(raw))

dat <- raw %>%
  dplyr::mutate(
    # --- outcome -----------------------------------------------------------
    asthma = asthma_from_items(MCQ010, MCQ035),

    # --- exposure ----------------------------------------------------------
    # SPISE = 600 * HDL-C^0.185 / (TG^0.2 * BMI^1.338), mg/dL and kg/m^2
    spise = 600 * LBDHDD^0.185 / (LBXTR^0.2 * BMXBMI^1.338),

    # --- covariates --------------------------------------------------------
    # "Don't know" (7) and "refused" (9) are set to missing for the
    # behavioural items, as in the published analysis.
    SMQ020 = dplyr::if_else(SMQ020 %in% c(7, 9), NA_real_, SMQ020),
    ALQ110 = dplyr::if_else(ALQ110 %in% c(7, 9), NA_real_, ALQ110),

    # Cardiovascular disease: any of the four self-reported conditions. A
    # "refused" (9) or "don't know" (7) answer on any item makes the composite
    # unclassifiable, so CVD is missing and the complete-case filter drops the
    # participant, as in the published analysis.
    CVD = dplyr::case_when(
      MCQ160B %in% c(7, 9) | MCQ160C %in% c(7, 9) |
        MCQ160D %in% c(7, 9) | MCQ160F %in% c(7, 9) ~ NA_integer_,
      MCQ160B == 1 | MCQ160C == 1 | MCQ160D == 1 | MCQ160F == 1 ~ 1L,
      TRUE ~ 0L
    ),
    cancer = dplyr::case_when(MCQ220 == 1 ~ 1L,
                              MCQ220 == 2 ~ 0L,
                              TRUE ~ NA_integer_),

    # Survey weight for the pooled five-cycle file
    WTMEC_combind = WTMEC2YR / WEIGHT_DIVISOR
  ) %>%
  # --- eligibility ---------------------------------------------------------
  dplyr::filter(RIDAGEYR >= AGE_MIN) %>%
  dplyr::filter(!is.na(asthma), !is.na(spise)) %>%
  # complete-case analysis across every covariate used in any model
  dplyr::filter(dplyr::if_all(
    c(RIAGENDR, RIDAGEYR, RIDRETH1, SMQ020, ALQ110, CVD, cancer,
      LBXGLU, LBXGH, LBDMONO, LBDLYMNO, LBDNENO, LBXRBCSI, LBXRDW, LBXPLTSI,
      LBDEONO, LBXTC, LBDLDL, WTMEC_combind),
    ~ !is.na(.x)
  )) %>%
  dplyr::mutate(
    # modelling copies with readable names
    female  = as.integer(RIAGENDR == 2),
    smoking = as.integer(SMQ020 == 1),
    drinking = as.integer(ALQ110 == 1),
    lymph   = LBDLYMNO,
    neut    = LBDNENO,
    eos     = LBDEONO,
    rdw     = LBXRDW,
    plt     = LBXPLTSI
  ) %>%
  dplyr::select(
    SEQN, cycle, asthma, spise, WTMEC_combind,
    RIAGENDR, RIDAGEYR, RIDRETH1, SMQ020, ALQ110, CVD, cancer,
    LBXGLU, LBXIN, LBXGH, LBDMONO, LBDLYMNO, LBDNENO, LBXRBCSI, LBXRDW,
    LBXPLTSI, LBDEONO,
    LBXTC, LBDLDL, LBDHDD, LBXTR, BMXBMI,
    SDMVPSU, SDMVSTRA,
    female, smoking, drinking, lymph, neut, eos, rdw, plt
  )

# ---------------------------------------------------------------------------
# 4. Checks
# ---------------------------------------------------------------------------
EXPECTED_N      <- 4588L
EXPECTED_EVENTS <- 413L

n <- nrow(dat); ev <- sum(dat$asthma)
message(sprintf("Analytic cohort: n = %d, events = %d (%.1f%%)", n, ev, 100 * ev / n))

if (n != EXPECTED_N || ev != EXPECTED_EVENTS) {
  warning(sprintf(
    paste0("Cohort size does not match the published analysis ",
           "(expected n = %d, events = %d; got n = %d, events = %d).\n",
           "Check that the raw files are the 2009-2018 releases and that the ",
           "eligibility rules in section 3 have not been edited."),
    EXPECTED_N, EXPECTED_EVENTS, n, ev), call. = FALSE)
} else {
  message("Cohort matches the published sample size.")
}

# ---------------------------------------------------------------------------
# 5. Development / validation split
# ---------------------------------------------------------------------------
# The split is by survey cycle, decided before any modelling: everything up to
# 2016 develops the model, 2017-2018 is held out for temporal validation.
split <- list(
  development = dplyr::filter(dat, cycle != "2017-2018"),
  validation  = dplyr::filter(dat, cycle == "2017-2018")
)

message(sprintf("Development (2009-2016): n = %d, events = %d",
                nrow(split$development), sum(split$development$asthma)))
message(sprintf("Validation  (2017-2018): n = %d, events = %d",
                nrow(split$validation), sum(split$validation$asthma)))

# ---------------------------------------------------------------------------
# 6. Save
# ---------------------------------------------------------------------------
write_step(dat, "cohort")
write_step(split, "split")

# A plain-text copy makes the derived dataset easy to inspect without R.
utils::write.csv(dat, path_data("cohort_analytic.csv"), row.names = FALSE)

stamp_session("02_build_cohort")
message("Wrote data/cohort.rds, data/split.rds, data/cohort_analytic.csv")
