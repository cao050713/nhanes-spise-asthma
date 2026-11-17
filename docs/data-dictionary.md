# Data dictionary

Every variable in `data/cohort_analytic.csv` (produced by
`R/02_build_cohort.R`), with its NHANES source. Variable names are the NHANES
item codes except where a derived name is used.

## Identifiers and design

| Variable | Source file | Description |
|---|---|---|
| `SEQN` | DEMO | Respondent sequence number |
| `cycle` | derived | Survey cycle, `2009-2010` … `2017-2018` |
| `SDMVPSU` | DEMO | Masked variance pseudo-PSU |
| `SDMVSTRA` | DEMO | Masked variance pseudo-stratum |
| `WTMEC2YR` | DEMO | Full-sample MEC exam weight, two-year |
| `WTMEC_combind` | derived | `WTMEC2YR / 5`, the pooled five-cycle weight |

## Outcome

| Variable | Source | Description |
|---|---|---|
| `MCQ010` | MCQ | "Ever told by a doctor you have asthma?" 1 yes, 2 no |
| `MCQ035` | MCQ | "Do you still have asthma?" 1 yes, 2 no (asked only when MCQ010 = 1) |
| `asthma` | derived | **Current asthma**: 1 when `MCQ010 = 1` and `MCQ035 = 1`; 0 when `MCQ010 = 2`, or `MCQ010 = 1` with `MCQ035 = 2`; otherwise missing. See `docs/KNOWN_DIFFERENCES.md` B1. |

## Exposure

| Variable | Source | Description |
|---|---|---|
| `LBDHDD` | HDL | HDL cholesterol, mg/dL |
| `LBXTR` | TRIGLY | Triglycerides, mg/dL |
| `BMXBMI` | BMX | Body mass index, kg/m² |
| `spise` | derived | `600 × LBDHDD^0.185 / (LBXTR^0.2 × BMXBMI^1.338)` |

## Demographic and behavioural covariates

| Variable | Source | Description |
|---|---|---|
| `RIAGENDR` | DEMO | 1 male, 2 female |
| `RIDAGEYR` | DEMO | Age in years at screening; 80 means 80 or older |
| `RIDRETH1` | DEMO | Race/ethnicity: 1 Mexican American, 2 Other Hispanic, 3 Non-Hispanic White, 4 Non-Hispanic Black, 5 Other |
| `SMQ020` | SMQ | Smoked at least 100 cigarettes in life: 1 yes, 2 no, 7 refused, 9 don't know |
| `ALQ110` | ALQ | Had at least 12 drinks in any one year: 1 yes, 2 no (item is `ALQ111` in 2017–2018) |
| `female` | derived | `RIAGENDR == 2` |
| `smoking` | derived | `SMQ020 == 1` |
| `drinking` | derived | `ALQ110 == 1` |

## Comorbidities

| Variable | Source | Description |
|---|---|---|
| `MCQ160B` | MCQ | Congestive heart failure |
| `MCQ160C` | MCQ | Coronary heart disease |
| `MCQ160D` | MCQ | Angina |
| `MCQ160F` | MCQ | Stroke |
| `CVD` | derived | 1 if any of the four above is 1; 0 if all are 2; missing if any is 7 or 9 |
| `MCQ220` | MCQ | Ever told you had cancer or malignancy: 1 yes, 2 no |
| `cancer` | derived | `MCQ220 == 1` |

## Laboratory covariates

| Variable | Source | Description | Unit |
|---|---|---|---|
| `LBXGLU` | GLU | Fasting glucose | mg/dL |
| `LBXGH` | GHB | Glycated haemoglobin | % |
| `LBDMONO` | CBC | Monocyte count | 1000 cells/µL |
| `LBDLYMNO` | CBC | Lymphocyte count | 1000 cells/µL |
| `LBDNENO` | CBC | Segmented neutrophil count | 1000 cells/µL |
| `LBDEONO` | CBC | Eosinophil count | 1000 cells/µL |
| `LBXRDW` | CBC | Red cell distribution width | % |
| `LBXPLTSI` | CBC | Platelet count | 1000 cells/µL |
| `LBXRBCSI` | CBC | Red blood cell count | million cells/µL |
| `LBXTC` | TCHOL | Total cholesterol | mg/dL |
| `LBDLDL` | TRIGLY | LDL cholesterol | mg/dL |
| `LBXIN` | INS | Fasting insulin (2013 onward only) | µIU/mL |

Short aliases used by the models: `lymph` = `LBDLYMNO`, `neut` = `LBDNENO`,
`eos` = `LBDEONO`, `rdw` = `LBXRDW`, `plt` = `LBXPLTSI`.

## Derived indices (created in `R/07_index_comparison.R`)

| Variable | Formula |
|---|---|
| `HOMA_IR` | `LBXIN × LBXGLU / 405` |
| `HOMA_BETA` | `360 × LBXIN / (LBXGLU − 63)`; participants with `LBXGLU ≤ 63` are excluded |
| `TyG` | `log(LBXTR × LBXGLU / 2)` |

## Eligibility rules applied in `R/02`

1. `RIDAGEYR >= 20`.
2. Non-missing `asthma` and non-missing `spise`.
3. Complete cases across every covariate used in any model: `RIAGENDR`,
   `RIDAGEYR`, `RIDRETH1`, `SMQ020`, `ALQ110`, `CVD`, `cancer`, `LBXGLU`,
   `LBXGH`, `LBDMONO`, `LBDLYMNO`, `LBDNENO`, `LBXRDW`, `LBXPLTSI`, `LBDEONO`,
   `LBXTC`, `LBDLDL`, `WTMEC_combind`.

Result: 4588 participants, 413 events.

## Coding conventions

* NHANES uses 7 for "refused" and 9 for "don't know"; both are set to missing
  for `SMQ020` and `ALQ110` before the complete-case filter.
* `RIDAGEYR` is top-coded at 80 in the public files. It is used as a
  continuous covariate.
* Survey weights are divided by the number of pooled cycles, as recommended for
  combined NHANES files.
