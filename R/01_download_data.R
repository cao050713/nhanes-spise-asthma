# =============================================================================
# 01_download_data.R -- fetch the NHANES files used in this study
#
# Downloads the public-use .XPT files for the five cycles of the National
# Health and Nutrition Examination Survey (2009-2010 through 2017-2018) into
# data/raw/. Files already present are skipped, so the script is safe to re-run.
#
#   Rscript R/01_download_data.R
#
# If the download fails (the CDC server sometimes throttles automated
# requests), fetch the files by hand from
#     https://wwwn.cdc.gov/nchs/nhanes/continuousnhanes/default.aspx?BeginYear=2017
# and drop them into data/raw/ with the same names.
# =============================================================================

source("R/utils.R")

# cycle suffix -> data release year used in the URL
cycles <- c(F = 2009, G = 2011, H = 2013, I = 2015, J = 2017)

# component -> variables we need from it (documentation only)
components <- list(
  DEMO   = "age, sex, race/ethnicity, survey design variables, MEC weight",
  BMX    = "body mass index",
  TRIGLY = "triglycerides (LBXTR), LDL cholesterol (LBDLDL)",
  HDL    = "HDL cholesterol (LBDHDD)",
  TCHOL  = "total cholesterol (LBXTC)",
  GLU    = "fasting glucose (LBXGLU)",
  GHB    = "glycated haemoglobin (LBXGH)",
  CBC    = "monocyte, lymphocyte, neutrophil, eosinophil, RDW, platelet, RBC",
  MCQ    = "asthma items (MCQ010, MCQ035), cancer (MCQ220), CVD items",
  SMQ    = "smoking (SMQ020)",
  ALQ    = "alcohol (ALQ110 in 2009-2016, ALQ111 in 2017-2018)",
  INS    = "fasting insulin (LBXIN); only available from 2013 onward"
)

base_url <- function(year, file, suffix) {
  sprintf("https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/%d/DataFiles/%s_%s.XPT",
          year, file, suffix)
}

# Insulin was not measured in 2009-2012.
components_cycle <- function(comp, suffix) {
  if (comp == "INS" && suffix %in% c("F", "G")) return(FALSE)
  TRUE
}

fetch <- function(url, dest, attempts = 3) {
  for (i in seq_len(attempts)) {
    ok <- tryCatch({
      utils::download.file(url, dest, mode = "wb", quiet = TRUE,
                           method = "libcurl")
      file.exists(dest) && file.size(dest) > 5000
    }, error = function(e) FALSE, warning = function(w) FALSE)

    if (isTRUE(ok)) return(TRUE)
    if (file.exists(dest)) file.remove(dest)
    Sys.sleep(2 * i)          # back off; the server throttles bursts
  }
  FALSE
}

manifest <- expand.grid(component = names(components), suffix = names(cycles),
                        stringsAsFactors = FALSE)
manifest$year <- cycles[manifest$suffix]
manifest <- manifest[!is.na(manifest$year), ]
manifest <- manifest[mapply(components_cycle, manifest$component, manifest$suffix), ]
manifest$file <- sprintf("%s_%s.XPT", manifest$component, manifest$suffix)
manifest$url  <- mapply(base_url, manifest$year, manifest$component,
                        manifest$suffix, USE.NAMES = FALSE)

message(sprintf("Need %d files (%d components x 5 cycles, minus insulin for 2009-2012).",
                nrow(manifest), length(components)))

failed <- character(0)
for (i in seq_len(nrow(manifest))) {
  dest <- path_raw(manifest$file[i])
  if (file.exists(dest) && file.size(dest) > 5000) {
    message(sprintf("[%2d/%2d] %-14s already present", i, nrow(manifest),
                    manifest$file[i]))
    next
  }
  message(sprintf("[%2d/%2d] %-14s downloading ...", i, nrow(manifest),
                  manifest$file[i]))
  if (!fetch(manifest$url[i], dest)) {
    failed <- c(failed, manifest$file[i])
    message("          FAILED")
  } else {
    message(sprintf("          ok (%.1f MB)", file.size(dest) / 1e6))
  }
}

utils::write.csv(manifest, path_data("file_manifest.csv"), row.names = FALSE)

if (length(failed)) {
  message("\n", length(failed), " file(s) could not be downloaded:\n  ",
          paste(failed, collapse = "\n  "))
  message("\nDownload them manually from the NHANES continuous survey page ",
          "for the matching cycle and place them in data/raw/:\n",
          "  https://wwwn.cdc.gov/nchs/nhanes/continuousnhanes/default.aspx")
  message("Then re-run this script; existing files are skipped.")
} else {
  message("\nAll files are in ", path_raw())
}
