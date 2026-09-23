# ==============================================================================
# 04  Geometric coupling between the response and the predictors
#
# Manuscript: Assessing the Landscape Composition and Configuration of
#             Threatened Medium- and Large-Sized Mammals using Habitat
#             Envelopes as Sampling Windows
# Author:     Maria Eduarda Nacif
# Study:      landscape configuration and the extent of the occupied habitat
#             envelope of nine threatened medium and large mammal taxa in the
#             Atlantic Forest of Rio de Janeiro State, Brazil.
#
# PURPOSE
#   Verify numerically, on all 67 units, the identity that ties patch density to
#   the response, log(extent) = log(number of patches) + log(100) - log(patch
#   density), and fit the reference model beside one legitimate control and two
#   circular ones. The circular controls are reported precisely so that the
#   circularity is visible rather than implied.
#
# INPUT   Dados/Processados/Data_Raw_WithCoords.csv
#         The same file that scripts 03 and 05 read. Until 23 August 2026
#         this script read Data_Raw_FINAL.csv, in which the response is
#         stored as an integer. The identity check below was therefore
#         computed on a rounded response and returned 0.0290, while script
#         05 returned 0.0277 for the same quantity on the full-precision
#         file. Reading the same input removes that discrepancy.
# OUTPUT  Outputs/Manuscrito/TableS28_Geometric_Coupling.csv
# ==============================================================================

suppressPackageStartupMessages({
  library(here)          # project-relative paths
  library(readr)         # read_csv
  library(dplyr)         # data manipulation
  library(lme4)          # linear mixed models
  library(rlang)         # sym(), used by the column-name fallback
})

GLOBAL_SEED <- 123
set.seed(GLOBAL_SEED)

data_path   <- here("Dados", "Processados")
output_path <- here("Outputs", "Manuscrito")
if (!dir.exists(output_path)) dir.create(output_path, recursive = TRUE)

message("\n[1/5] Loading data...")

# Data_Raw_WithCoords.csv is written by script 01 and carries the response at
# full precision, together with the repaired class metrics. Data_Raw_FINAL.csv is
# an earlier file kept only so that an old run can still be reproduced; it stores
# the response as an integer and must not be used here.
primary <- file.path(data_path, "Data_Raw_WithCoords.csv")
if (!file.exists(primary)) {
  stop("Data_Raw_WithCoords.csv not found in ", data_path,
       ". Run scripts 00 and 01 before this one.", call. = FALSE)
}
d <- read_csv(primary, show_col_types = FALSE)

# Back-compatibility with the AOH_ha naming used in Part III. AOH_unit_area_ha is
# the name script 01 writes; the other two are earlier names.
for (nm in c("AOH_unit_area_ha", "UA_Area_ha", "EHA_ha")) {
  if (nm %in% names(d) && !"AOH_ha" %in% names(d)) {
    d <- d %>% rename(AOH_ha = !!rlang::sym(nm))
  }
}
if (!"AOH_ha" %in% names(d)) {
  stop("No response column found in Data_Raw_WithCoords.csv. Expected one of ",
       "AOH_ha, AOH_unit_area_ha, UA_Area_ha or EHA_ha.", call. = FALSE)
}

d <- d %>%
  mutate(log_AOH = log(AOH_ha),
         log_NP  = log(NP_Forest),
         log_CA  = log(pmax(CA_Forest, 1)),
         GENUS   = as.factor(GENUS))

# Guard: on the full-precision response the corrected dataset has
# mean log(AOH) = 10.3735 and sd = 0.9392. On the rounded response it had
# 10.3758 and 0.9508, and the version carrying the join error had 10.3862.
# A tolerance of 0.01 accepts the first two and rejects the third.
message(sprintf("   n = %d | mean log(OHE) = %.4f | sd = %.4f",
                nrow(d), mean(d$log_AOH), sd(d$log_AOH)))
if (abs(mean(d$log_AOH) - 10.3735) > 0.01) {
  warning("mean log(OHE) does not match the corrected dataset. ",
          "Confirm that script 01 has been run and that ",
          "Data_Raw_WithCoords.csv is the file it wrote, before ",
          "interpreting this table.")
}

# Standardize every variable that enters a model, so coefficients are comparable.
z_vars <- c("Mean_FRAC_Forest", "PD_Forest", "PROX_MN_Forest",
            "PLAND_Forest", "log_NP", "log_CA")
d <- d %>% mutate(across(all_of(z_vars), ~as.numeric(scale(.x)), .names = "{.col}_z"))

# ------------------------------------------------------------------------------
# 2. The four models
# ------------------------------------------------------------------------------
message("[2/5] Fitting the four mixed models...")

base3 <- c("Mean_FRAC_Forest_z", "PD_Forest_z", "PROX_MN_Forest_z")

fit_mixed <- function(preds) {
  f <- as.formula(paste("log_AOH ~", paste(preds, collapse = " + "), "+ (1 | GENUS)"))
  # REML = FALSE so that the models remain comparable by likelihood.
  lme4::lmer(f, data = d, REML = FALSE)
}

model_set <- list(
  list(name = "M1. Retained model",
       preds = base3,
       note  = "reference"),
  list(name = "M2. + PLAND (proportion of forest)",
       preds = c(base3, "PLAND_Forest_z"),
       note  = "legitimate control: PLAND is scale-free"),
  list(name = "M3. + log(NP) (number of patches)",
       preds = c(base3, "log_NP_z"),
       note  = "circular: log(OHE) = log(NP) + log(100) - log(PD)"),
  list(name = "M4. + log(CA) (absolute forest area)",
       preds = c(base3, "log_CA_z"),
       note  = "circular: CA = PLAND x AOH / 100")
)

rows <- list()

for (spec in model_set) {
  m  <- fit_mixed(spec$preds)
  cf <- summary(m)$coefficients          # Estimate, Std. Error, t value
  # lme4 does not return p-values by design; the Wald normal approximation is
  # used here, which is the same convention adopted in Part III.
  for (v in spec$preds) {
    est <- cf[v, 1]; se <- cf[v, 2]; tv <- cf[v, 3]
    rows[[length(rows) + 1]] <- data.frame(
      Model     = spec$name,
      Predictor = sub("_z$", "", v),
      Estimate  = round(est, 4),
      SE        = round(se, 4),
      P_value   = round(2 * pnorm(abs(tv), lower.tail = FALSE), 4),
      Note      = spec$note,
      stringsAsFactors = FALSE
    )
  }
  vc      <- as.data.frame(lme4::VarCorr(m))
  var_genus   <- vc$vcov[vc$grp == "GENUS"]
  var_resid   <- vc$vcov[vc$grp == "Residual"]
  rows[[length(rows) + 1]] <- data.frame(
    Model = spec$name, Predictor = "ICC (genus)",
    Estimate = round(var_genus / (var_genus + var_resid), 4),
    SE = NA_real_, P_value = NA_real_, Note = spec$note,
    stringsAsFactors = FALSE
  )
  message(sprintf("   %-38s ICC = %.4f", spec$name, var_genus / (var_genus + var_resid)))
}

# ------------------------------------------------------------------------------
# 3. Partial correlations of PROX_MN with log(AOH)
# ------------------------------------------------------------------------------
message("[3/5] Partial correlations...")

partial_cor <- function(x, y, ctrl) {
  if (length(ctrl) == 0) return(cor(d[[x]], d[[y]]))
  fx <- as.formula(paste(x, "~", paste(ctrl, collapse = " + ")))
  fy <- as.formula(paste(y, "~", paste(ctrl, collapse = " + ")))
  cor(residuals(lm(fx, data = d)), residuals(lm(fy, data = d)))
}

controls <- list(
  list(v = character(0),                      lab = "none",            note = "legitimate control"),
  list(v = "PLAND_Forest",                    lab = "PLAND",           note = "legitimate control"),
  list(v = "log_NP",                          lab = "log(NP)",         note = "circular control"),
  list(v = c("PLAND_Forest", "log_NP"),       lab = "PLAND + log(NP)", note = "circular control")
)

for (ctl in controls) {
  r <- partial_cor("PROX_MN_Forest", "log_AOH", ctl$v)
  rows[[length(rows) + 1]] <- data.frame(
    Model = "Partial correlation of PROX_MN with log(OHE)",
    Predictor = paste("controlling for", ctl$lab),
    Estimate = round(r, 4), SE = NA_real_, P_value = NA_real_,
    Note = ctl$note, stringsAsFactors = FALSE
  )
  message(sprintf("   controlling for %-16s r = %.4f", ctl$lab, r))
}

# ------------------------------------------------------------------------------
# 4. Numerical check of the algebraic identity
# ------------------------------------------------------------------------------
message("[4/5] Algebraic identity check...")

max_discrepancy <- max(abs(log(d$NP_Forest) + log(100) - log(d$PD_Forest) - d$log_AOH))
rows[[length(rows) + 1]] <- data.frame(
  Model = "Algebraic identity check",
  Predictor = "max |log(NP)+log(100)-log(PD) - log(OHE)|",
  Estimate = round(max_discrepancy, 4), SE = NA_real_, P_value = NA_real_,
  Note = "exact identity: PD is defined as NP/(OHE/100)",
  stringsAsFactors = FALSE
)
message(sprintf("   maximum absolute discrepancy across %d units: %.4f", nrow(d), max_discrepancy))

# ------------------------------------------------------------------------------
# 5. Export
# ------------------------------------------------------------------------------
message("[5/5] Exporting...")

out_table <- do.call(rbind, rows)
write.csv(out_table, file.path(output_path, "TableS28_Geometric_Coupling.csv"),
          row.names = FALSE, na = "")

message("\n============================================================")
message("  PART IV - Geometric coupling: summary")
message("============================================================")
print(out_table, row.names = FALSE)
message("\n  Exported: TableS28_Geometric_Coupling.csv")
message("  Read M2 as the test. M3 and M4 are circular by construction and are")
message("  reported only so that the circularity is visible to the reader.")
message("============================================================\n")

message("\n[Session]")
print(sessionInfo())
writeLines(capture.output(sessionInfo()),
           file.path(output_path, "session_info_part04.txt"))

