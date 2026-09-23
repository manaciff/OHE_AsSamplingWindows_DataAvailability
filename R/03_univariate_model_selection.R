# ==============================================================================
# 03  Univariate model of envelope extent: specification, inference and figures
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
#   Establish which forest configuration metrics accompany a larger envelope,
#   under a candidate set defined a priori.
#
#   The goal is inference, not prediction. Six specifications are compared on a
#   common response scale by AICc, by BIC and by exact leave-one-out cross
#   validation, with leave-one-genus-out beside it. Both information criteria
#   are reported because their relative behavior depends on unobserved
#   heterogeneity, which 67 units from nine taxa carry (Brewer et al., 2016).
#   The full a priori model is then fitted once, without selection, and its
#   coefficients are the inferential result: inference on coefficients is biased
#   when it follows a selection step (Yates et al., 2023; Tredennick et al.,
#   2021). The same model is refitted with the Moran eigenvectors retained by
#   script 02, because an information criterion computed on a non-spatial
#   regression is sensitive to spatial autocorrelation (Diniz-Filho et al.,
#   2008). The exhaustive search over subsets is reported as a description of
#   selection uncertainty, not as a second set of estimates to test.
#
#   The interaction between proximity and forest cover is the only one admitted,
#   because the fragmentation threshold hypothesis is a statement about that
#   interaction rather than about either term alone.
#
# INPUT   Dados/Processados/Data_Raw_WithCoords.csv; TableII_2b_MEM_Vectors.csv
# OUTPUT  Outputs/Manuscrito/Table_S_PartIII_*.csv and Figure_PIII_*.tiff.
#         Table_S_PartIII_ConfirmatoryModel.csv carries the inferential
#         estimates; Table_S_PartIII_SpatialSensitivity.csv shows whether they
#         hold with spatial structure in the model; the averaging table
#         describes selection uncertainty only.
#
# IMPORTANT
#   The model selected here is not the model the manuscript interprets. Script
#   05 explains why. This script is retained because the candidate set was
#   pre-specified and reporting the pre-specified analysis is part of the record.
# ==============================================================================


# ==============================================================================
# 1. ENVIRONMENT SETUP
# ==============================================================================
suppressPackageStartupMessages({
  library(here)        # reproducible paths
  library(tidyverse)   # data wrangling and plotting
  library(scales)      # axis formatting
  library(patchwork)   # multi-panel composition
  library(MuMIn)       # multimodel inference (dredge, model.avg)
  library(DHARMa)      # simulation-based residual diagnostics
  library(car)         # vif, ncvTest (Breusch-Pagan), durbinWatsonTest
  library(lme4)        # linear mixed models (lmer), one of the six candidates
  library(performance) # icc(), r2() for mixed models (variance components)
})

# Reproducibility seed; replicated locally before each stochastic call.
GLOBAL_SEED <- 123

# Minimum number of observations a genus must have for a per-genus association
# to be estimated in the heatmap of Section 10.5. Below this, the cell is shown
# as "not estimable" (grey) rather than reporting an unstable correlation.
N_UNITS_EXPECTED <- 67   # the design has 67 sampling units; Section 3 asserts it
MIN_N_GENUS <- 4

# Colourblind-safe diverging palette for the association heatmap (Section 10.5).
# Blue for negative associations, near-white at zero, vermillion for positive,
# consistent with the Okabe-Ito family used across Parts I-III.
heat_low  <- "#0072B2"   # blue       (negative association)
heat_mid  <- "#F7F7F7"   # near-white (no association)
heat_high <- "#D55E00"   # vermillion (positive association)

# Configurable paths, relative to the project root.
data_path   <- here("Dados", "Processados")
output_path <- here("Outputs", "Manuscrito")

# Ensure output directory exists.
if (!dir.exists(output_path)) dir.create(output_path, recursive = TRUE)

# ------------------------------------------------------------------------------
# 1.1 Color palette for Part III figures (Springer-friendly, colorblind safe)
# ------------------------------------------------------------------------------
# A single Okabe-Ito-derived palette is used across all Part III figures to
# keep visual identity consistent with Parts I and II.
colors_predictors <- c(
  "PD_Forest"            = "#009E73",   # green       (patch density)
  "Mean_FRAC_Forest"     = "#0072B2",   # blue        (patch shape)
  "Mean_ENN_Forest_m"    = "#CC79A7",   # rose        (isolation)
  "PROX_MN_Forest"       = "#E69F00"    # orange      (forest proximity)
)

# Human-readable predictor labels used in figures and tables. Each label
# mirrors the FRAGSTATS/landscapemetrics nomenclature so that the figures can
# be read without consulting the methods section.
predictor_labels <- c(
  "PD_Forest"            = "Forest patch\ndensity (PD)",
  "Mean_FRAC_Forest"     = "Forest patch\nshape (mean FRAC)",
  "Mean_ENN_Forest_m"    = "Forest patch\nisolation (mean ENN)",
  "PROX_MN_Forest"       = "Forest patch\nproximity (PROX_MN)",
  "PLAND_Forest"         = "Forest cover\n(PLAND)",
  "PLAND_Forest:PROX_MN_Forest" = "Forest cover x\nproximity",
  "PROX_MN_Forest:PLAND_Forest" = "Forest cover x\nproximity"
)

# The standardization suffix has to come off every term of an interaction, not
# only off the last one: "PLAND_Forest_z:PROX_MN_Forest_z" must reduce to
# "PLAND_Forest:PROX_MN_Forest". Anchoring on "_z$" alone leaves the interior
# suffix in place and the lookup then misses.
strip_z <- function(x) gsub("_z(?=:|$)", "", x, perl = TRUE)

# Look up a label and fall back to the raw variable name when the lookup misses.
# A bare `predictor_labels[x]` returns NA for any name absent from the vector,
# and an NA propagates silently into a factor level, which is how an axis label
# reading "NA" reached a figure. This wrapper makes a missing entry visible as
# the variable name instead of as a hole, and warns once so the omission is
# fixed at the source rather than tolerated.
label_predictors <- function(x) {
  out <- unname(predictor_labels[x])
  gap <- is.na(out)
  if (any(gap)) {
    warning("No label defined for: ", paste(unique(x[gap]), collapse = ", "),
            ". Falling back to the variable name.", call. = FALSE)
    out[gap] <- x[gap]
  }
  out
}

# ==============================================================================
# 2. ANALYTICAL UTILITIES
# ==============================================================================

# Scalar-safe null/NA-coalescing operator. Returns `b` if `a` is NULL, of length
# Safe scalar coercion: extracts the first finite numeric element from `x`,
# returning NA_real_ when `x` is NULL, empty, or non-numeric. Used to tame the
# heterogeneous return shapes of car:: and DHARMa:: test objects before
# rounding.
safe_num <- function(x) {
  if (is.null(x) || length(x) == 0) return(NA_real_)
  v <- suppressWarnings(as.numeric(x[[1]]))
  if (!is.finite(v)) return(NA_real_)
  v
}

# A) Publication-quality plot export (TIFF, 600 dpi, LZW compression)
# Reports the absolute path and modification time of every figure it writes, so
# it is always clear WHERE and WHEN each figure was (re)generated. A single
# failed write (most often a file left open in an image viewer, which locks it
# against overwrite) is reported as a warning and does not halt the run.
save_publication_plot <- function(plot_obj, file_name,
                                  width_mm = 174, height_mm = 130) {
  if (is.null(plot_obj)) {
    message(sprintf("   [skip] %s: plot object is NULL.", file_name))
    return(invisible(FALSE))
  }
  full_path <- file.path(output_path, file_name)
  ok <- tryCatch({
    ggsave(filename = file_name, plot = plot_obj, path = output_path,
           width = width_mm, height = height_mm, units = "mm",
           dpi = 600, device = "tiff", compression = "lzw", bg = "white")
    TRUE
  }, error = function(e) {
    message(sprintf(paste0(
      "   [WARN] Could not write %s\n",
      "          Reason: %s\n",
      "          If the file is open in an image viewer (Photos, Preview, etc.), ",
      "close it and re-run: a locked file cannot be overwritten."),
      normalizePath(full_path, mustWork = FALSE), conditionMessage(e)))
    FALSE
  })
  if (isTRUE(ok) && file.exists(full_path)) {
    message(sprintf("   Figure written: %s  (%s)",
                    normalizePath(full_path, mustWork = FALSE),
                    format(file.info(full_path)$mtime, "%Y-%m-%d %H:%M:%S")))
  }
  invisible(ok)
}

# B) Compact residual diagnostic battery for OLS and GLM
# Returns a tidy data frame with one row per test, plus a flag indicating
# whether all assumptions are met at alpha = 0.05.
run_residual_diagnostics <- function(model, label = "model", seed = GLOBAL_SEED) {

  # WHICH TESTS APPLY TO WHICH FAMILY (added 6 Aug 2026)
  # Shapiro-Wilk, the score test of non-constant variance and Durbin-Watson are
  # tests of the Gaussian assumptions: normal residuals with constant variance.
  # A Gamma GLM makes neither assumption, so applying them to its raw residuals
  # and reporting "assumption not met" is a category error: it flags the model
  # for failing something it never claimed. Those three are therefore reported
  # as not applicable outside the Gaussian families, and the simulation-based
  # diagnostics, which are defined for any family, carry the verdict.
  gaussiano <- inherits(model, "lm") && !inherits(model, "glm") ||
               inherits(model, "lmerMod")

  # Shapiro-Wilk on raw residuals.
  sw <- if (gaussiano) {
    tryCatch(shapiro.test(residuals(model)),
             error = function(e) list(statistic = NA_real_, p.value = NA_real_))
  } else list(statistic = NA_real_, p.value = NA_real_)
  
  # Breusch-Pagan / non-constant variance via car::ncvTest.
  bp <- if (gaussiano) {
    tryCatch(car::ncvTest(model),
             error = function(e) list(ChiSquare = NA_real_, p = NA_real_))
  } else list(ChiSquare = NA_real_, p = NA_real_)

  # Durbin-Watson for residual serial dependence.
  dw <- if (gaussiano) {
    tryCatch(car::durbinWatsonTest(model, reps = 999),
             error = function(e) list(dw = NA_real_, p = NA_real_))
  } else list(dw = NA_real_, p = NA_real_)
  
  # DHARMa simulation-based diagnostics.
  set.seed(seed)
  sim  <- DHARMa::simulateResiduals(fittedModel = model, n = 1000, seed = seed)
  disp <- DHARMa::testDispersion(sim, plot = FALSE)
  unif <- DHARMa::testUniformity(sim, plot = FALSE)
  outl <- DHARMa::testOutliers(sim, plot = FALSE)
  
  diag_df <- data.frame(
    Model     = label,
    Test      = c("Shapiro-Wilk (normality)",
                  "Breusch-Pagan (homoscedasticity)",
                  "Durbin-Watson (independence)",
                  "DHARMa dispersion",
                  "DHARMa uniformity",
                  "DHARMa outliers"),
    Statistic = round(c(safe_num(sw$statistic),
                        safe_num(bp$ChiSquare),
                        safe_num(dw$dw),
                        safe_num(disp$statistic),
                        safe_num(unif$statistic),
                        safe_num(outl$statistic)), 4),
    P_value   = round(c(safe_num(sw$p.value),
                        safe_num(bp$p),
                        safe_num(dw$p),
                        safe_num(disp$p.value),
                        safe_num(unif$p.value),
                        safe_num(outl$p.value)), 4),
    Assumption_met = NA,
    stringsAsFactors = FALSE
  )
  
  # A test that does not apply to the family is reported as such, never as a
  # failed assumption.
  diag_df$Assumption_met <- ifelse(is.na(diag_df$P_value), "not applicable",
                                   ifelse(diag_df$P_value > 0.05, "TRUE", "FALSE"))
  diag_df
}

# C) Coefficient extraction with 95% confidence intervals
# Works for lm, glm, and model.avg objects.
extract_coefs <- function(model, label = "final") {
  if (inherits(model, "averaging")) {
    # Conditional (subset) averages: parameter averaged only across models
    # where it occurs. Reported as default in Burnham & Anderson (2002).
    cm <- summary(model)$coefmat.subset
    ci <- confint(model, full = FALSE)
    df <- data.frame(
      Predictor = rownames(cm),
      Estimate  = cm[, "Estimate"],
      SE        = cm[, "Std. Error"],
      CI_lower  = ci[, 1],
      CI_upper  = ci[, 2],
      Z_value   = cm[, "z value"],
      P_value   = cm[, "Pr(>|z|)"],
      Source    = label,
      stringsAsFactors = FALSE
    )
  } else if (inherits(model, "merMod")) {
    # ADDED 29 Jul 2026. lme4 returns only three columns (Estimate, Std. Error,
    # t value) and no p-value, and confint() on a merMod prepends rows for the
    # random-effect standard deviations, so the generic branch below fails with
    # a subscript error. This is where the previous run stopped, immediately
    # after Table_S_PartIII_VariableImportance.csv was written. Fixed effects
    # only are extracted, and the p-value is the Wald normal approximation,
    # which is the convention consistent with the Z_value column used
    # throughout this script.
    coefs <- summary(model)$coefficients
    ci    <- suppressMessages(
      confint(model, parm = "beta_", method = "Wald", oldNames = FALSE))
    ci    <- ci[rownames(coefs), , drop = FALSE]
    tval  <- coefs[, 3]
    df <- data.frame(
      Predictor = rownames(coefs),
      Estimate  = coefs[, 1],
      SE        = coefs[, 2],
      CI_lower  = ci[, 1],
      CI_upper  = ci[, 2],
      Z_value   = tval,
      P_value   = 2 * stats::pnorm(abs(tval), lower.tail = FALSE),
      Source    = label,
      stringsAsFactors = FALSE
    )
  } else {
    coefs <- summary(model)$coefficients
    ci    <- suppressMessages(confint(model))
    df <- data.frame(
      Predictor = rownames(coefs),
      Estimate  = coefs[, 1],
      SE        = coefs[, 2],
      CI_lower  = ci[, 1],
      CI_upper  = ci[, 2],
      Z_value   = coefs[, 3],
      P_value   = coefs[, 4],
      Source    = label,
      stringsAsFactors = FALSE
    )
  }
  rownames(df) <- NULL
  df
}

# D) Deviance-based pseudo-R^2 for GLMs (Cohen et al. 2003)
# For lm, the adjusted R^2 from summary() is returned instead.
explained_variance <- function(model) {
  if (inherits(model, "merMod")) {
    # ADDED 30 Jul 2026. Marginal R2 of Nakagawa & Schielzeth (2013): the share
    # explained by the fixed effects alone, which is the quantity comparable to
    # the adjusted R2 of the OLS and the one reported in the manuscript.
    r2 <- tryCatch(performance::r2(model)$R2_marginal,
                   error = function(e) NA_real_)
    as.numeric(r2)
  } else if (inherits(model, "glm")) {
    1 - (model$deviance / model$null.deviance)
  } else if (inherits(model, "lm")) {
    summary(model)$adj.r.squared
  } else {
    NA_real_
  }
}

# ==============================================================================
# 3. DATA LOADING AND PREPARATION
# ==============================================================================
message("\n[Initialization] Loading dataset and preparing Part III predictors...")

# Predictor set retained from Part I (RDA term test), restricted to the
# Forest class for Part III (see header revision note v3-i).
predictors_part3 <- c("PD_Forest", "Mean_FRAC_Forest",
                      "Mean_ENN_Forest_m", "PROX_MN_Forest",
                      "PLAND_Forest")

# ------------------------------------------------------------------------------
# A PRIORI INTERACTION (added 6 Aug 2026)
# The fragmentation threshold hypothesis states that the spatial arrangement of
# habitat matters where habitat is scarce and ceases to matter where it is
# abundant (Andren 1994; Pardini et al. 2010; Villard & Metzger 2014). That is a
# statement about an interaction between connectivity and habitat amount, not
# about either term on its own, so the product of the two is admitted to the
# candidate set. It is the ONLY interaction considered: with 67 units, a full
# factorial of the four configuration metrics would leave fewer than eight
# observations per parameter.
#
# Predictors are z-standardized before the product is formed, which keeps the
# constituent coefficients interpretative at the mean and removes the
# non-essential collinearity that products introduce (Aiken & West 1991;
# Schielzeth 2010). dredge() is instructed to respect marginal, so no model
# containing the product is fitted without both constituent terms.
# --- A) Primary data file ---
# CORRECTED 22 Aug 2026. This block used to accept three files in order of
# preference: Data_Raw_WithCoords.csv, then Data_Raw_FINAL.csv, then
# Data_Raw_GLM_NMDS_Final.csv. Only the first carries the metric repair of
# Part 0 and the canonical taxon labels, so only the first reproduces the
# published results. The fallbacks were kept for back-compatibility and made
# the script silently reproducible to a different answer, which is the opposite
# of what a fallback should do. The single required input is now named, and the
# run stops with an instruction if it is absent.
# Same definition as Parts 0, I and II. Keep the four copies in step.
canonical_taxon <- function(x) {
  key <- gsub("[^a-z]", "", tolower(as.character(x)))
  out <- dplyr::case_when(
    key %in% c("aguariba", "aguriba", "alouattaguariba")           ~ "A_guariba",
    key %in% c("barach", "baracnoides", "barachnoides",
               "brachytelesarachnoides")                           ~ "B_arachnoides",
    key %in% c("bcrinitus", "btorquatus", "bradypuscrinitus")       ~ "B_crinitus",
    key %in% c("bhypox", "bhypoxantus", "bhypoxanthus",
               "brachyteleshypoxanthus")                           ~ "B_hypoxanthus",
    key %in% c("lwiedii", "leoparduswiedii")                        ~ "L_wiedii",
    key %in% c("mazama", "mazamasp")                                ~ "Mazama",
    key %in% c("mtrida", "mtridactyla", "myrmecophagatridactyla")   ~ "M_tridactyla",
    key %in% c("pconcolor", "concolor", "pumaconcolor")             ~ "P_concolor",
    key %in% c("tpecari", "tayassupecari")                          ~ "T_pecari",
    TRUE ~ NA_character_
  )
  if (any(is.na(out) & !is.na(x))) {
    stop("Unrecognised taxon name(s): ",
         paste(unique(x[is.na(out) & !is.na(x)]), collapse = ", "))
  }
  out
}

REQUIRED_INPUT <- "Data_Raw_WithCoords.csv"
primary_file   <- file.path(data_path, REQUIRED_INPUT)
if (!file.exists(primary_file)) {
  stop(sprintf(paste0("%s not found in %s. It is written by script 01 and is ",
                      "the only input that carries the metric repair of script ",
                      "00 and the canonical taxon labels. Run scripts 00 and 01 ",
                      "before this one. Data_Raw_FINAL.csv and ",
                      "Data_Raw_GLM_NMDS_Final.csv are earlier files and are no ",
                      "longer accepted here: they yield a different sample."),
               REQUIRED_INPUT, data_path))
}
message(sprintf("   Primary data file: %s", basename(primary_file)))

data_raw <- read_csv(primary_file, show_col_types = FALSE)

# Back-compatibility: rename legacy column names to the current AOH_ha standard.
# Both UA_Area_ha (legacy) and EHA_ha (a transitional name used in v1.x) are
# accepted as input and renamed to AOH_ha for internal consistency.
# AOH_unit_area_ha is the name Part I writes into Data_Raw_WithCoords.csv, which
# is now the primary input. It must be listed here, otherwise the response
# variable is not found and the run fails at the first model.
if ("AOH_unit_area_ha" %in% names(data_raw) && !"AOH_ha" %in% names(data_raw)) {
  data_raw <- data_raw %>% rename(AOH_ha = AOH_unit_area_ha)
}
if ("UA_Area_ha" %in% names(data_raw) && !"AOH_ha" %in% names(data_raw)) {
  data_raw <- data_raw %>% rename(AOH_ha = UA_Area_ha)
}
if ("EHA_ha" %in% names(data_raw) && !"AOH_ha" %in% names(data_raw)) {
  data_raw <- data_raw %>% rename(AOH_ha = EHA_ha)
}
if (!"AOH_ha" %in% names(data_raw)) {
  stop("No response column found. Expected one of AOH_ha, AOH_unit_area_ha, ",
       "UA_Area_ha or EHA_ha in ", basename(primary_file), ".")
}

# Canonical taxon labels, so that SPECIES and GENUS match Parts 0, I and II.
if ("SPECIES" %in% names(data_raw)) {
  data_raw$SPECIES <- canonical_taxon(data_raw$SPECIES)
}

# --- B) Auxiliary file: only bind columns that are GENUINELY ABSENT ---
# Avoids dplyr's suffix collision (.x / .y) when the primary file already
# contains the predictor. The legacy file is consulted only as a complement.
legacy_candidates <- c("UA_Variables_GLM.csv", "UA_Variables_GLM.CSV")
legacy_path <- NULL
for (cand in legacy_candidates) {
  candidate_path <- file.path(data_path, cand)
  if (file.exists(candidate_path)) {
    legacy_path <- candidate_path
    break
  }
}

missing_now <- setdiff(predictors_part3, names(data_raw))
if (length(missing_now) > 0 && !is.null(legacy_path)) {
  message(sprintf("   Auxiliary file detected: %s", basename(legacy_path)))
  message(sprintf("   Columns to be supplemented from auxiliary file: %s",
                  paste(missing_now, collapse = ", ")))
  
  # CORRECTED 29 Jul 2026. The previous version selected UA_ID as the only join
  # key and applied distinct(UA_ID, .keep_all = TRUE). UA_ID is a WITHIN-TAXON
  # polygon index: every taxon has its own sequence starting at 1. Joining on
  # UA_ID alone therefore propagated the values of the first taxon holding a
  # given index to every other taxon sharing it, which corrupted 52 cells across
  # 8 rows (elevation, human density and all PROX_MN columns of A_guariba
  # UA 10-12, Mazama UA 10 and P_concolor UA 10-12). The key must be the pair
  # (SPECIES, UA_ID), which identifies a sampling unit uniquely.
  dados_antigos <- read_csv(legacy_path, show_col_types = FALSE) %>%
    select(any_of(c("SPECIES", "UA_ID", missing_now)))

  if (!"SPECIES" %in% names(dados_antigos)) {
    stop("The auxiliary file must carry a SPECIES column; joining on UA_ID ",
         "alone mixes taxa and silently corrupts the predictors.")
  }
  if (anyDuplicated(dados_antigos[, c("SPECIES", "UA_ID")]) > 0) {
    stop("The auxiliary file has duplicated (SPECIES, UA_ID) pairs; resolve ",
         "them before joining.")
  }

  # Cast join keys to a common type to avoid silent join failures.
  dados_antigos$UA_ID   <- as.character(dados_antigos$UA_ID)
  dados_antigos$SPECIES <- as.character(dados_antigos$SPECIES)
  data_raw$UA_ID        <- as.character(data_raw$UA_ID)
  data_raw$SPECIES      <- as.character(data_raw$SPECIES)

  n_before <- nrow(data_raw)
  data_raw <- data_raw %>%
    left_join(dados_antigos, by = c("SPECIES", "UA_ID"))
  stopifnot(nrow(data_raw) == n_before)
} else if (length(missing_now) == 0) {
  message("   All Part III predictors already present in the primary file. ",
          "Auxiliary file not required.")
}

# Coerce categorical columns to factor.
data_raw <- data_raw %>%
  mutate(across(any_of(c("ORDER", "GENUS", "SPECIES",
                         "DIET", "LOCOMOTION", "UA_ID")), as.factor))

# Defensive check: stop with a clear message if any predictor is still absent.
missing_preds <- setdiff(predictors_part3, names(data_raw))
if (length(missing_preds) > 0) {
  stop(sprintf(paste0("The following predictor(s) are missing from the input ",
                      "data after merging: %s. Confirm that Data_Raw_FINAL.csv ",
                      "(or the auxiliary UA_Variables_GLM file) contains these ",
                      "columns before running Part III."),
               paste(missing_preds, collapse = ", ")))
}

# Working data frame for Part III. A small constant, 0.001 ha, is added to any
# observation with AOH = 0 so that the logarithm is defined. On these data no
# unit has an area of zero, so the constant is never applied; it is kept as a
# guard. The same literal appears at the second occurrence of this operation,
# in Section 7.5.
data_p3 <- data_raw %>%
  select(UA_ID, GENUS, SPECIES, AOH_ha, all_of(predictors_part3)) %>%
  drop_na(all_of(predictors_part3)) %>%
  mutate(AOH_ha   = ifelse(AOH_ha == 0, 0.001, AOH_ha),
         log_AOH  = log(AOH_ha))

# Z-score predictors so that beta coefficients are directly comparable.
data_p3 <- data_p3 %>%
  mutate(across(all_of(predictors_part3),
                ~as.numeric(scale(.x)),
                .names = "{.col}_z"))

n_obs <- nrow(data_p3)
n_su  <- length(unique(data_p3$UA_ID))
message(sprintf("   Sample size after pre-processing: n = %d observations distributed across %d sampling units (SUs).",
                n_obs, n_su))

# ------------------------------------------------------------------------------
# GUARD ON THE SAMPLE SIZE (v7). drop_na() above is a guard, not a working step.
# The full design is 67 units; anything smaller means the input was not produced
# by Part I (v10 or later), or that Part 0 was not run. The run stops and names
# the offending cells rather than analysing a subset and reporting it in a
# supplementary table.
# ------------------------------------------------------------------------------
if (n_obs != N_UNITS_EXPECTED) {
  lost <- data_raw %>%
    select(SPECIES, UA_ID, all_of(predictors_part3)) %>%
    filter(!complete.cases(.)) %>%
    pivot_longer(-c(SPECIES, UA_ID), names_to = "Column", values_to = "Value") %>%
    filter(is.na(Value))
  if (nrow(lost) > 0) print(as.data.frame(lost), row.names = FALSE)
  stop(sprintf(paste0("Part III is running on %d observations instead of %d. ",
                      "Run Script_Part_0_Recompute_Units_v2.R, then Part I ",
                      "(v10 or later), and make sure Data_Raw_WithCoords.csv is ",
                      "the file being read."),
               n_obs, N_UNITS_EXPECTED))
}
message(sprintf("   Average of %.2f observations per SU.", n_obs / n_su))

# Export descriptive statistics for reporting in the Methods section.
desc_stats <- data_p3 %>%
  select(AOH_ha, log_AOH, all_of(predictors_part3)) %>%
  pivot_longer(everything(), names_to = "Variable", values_to = "Value") %>%
  group_by(Variable) %>%
  summarise(N      = sum(!is.na(Value)),
            Mean   = round(mean(Value, na.rm = TRUE), 4),
            SD     = round(sd(Value, na.rm = TRUE),   4),
            Median = round(median(Value, na.rm = TRUE), 4),
            Min    = round(min(Value, na.rm = TRUE),   4),
            Max    = round(max(Value, na.rm = TRUE),   4),
            .groups = "drop") %>%
  # The manuscript withdrew the term Area of Habitat, so the exported labels use
  # OHE. Only the printed labels change; the internal column names are untouched.
  mutate(Variable = ifelse(Variable == "AOH_ha",  "OHE_ha",
                    ifelse(Variable == "log_AOH", "log_OHE", Variable)))
write.csv(desc_stats,
          file.path(output_path, "Table_S_PartIII_DescriptiveStatistics.csv"),
          row.names = FALSE)
message("   Descriptive statistics exported.")

# ==============================================================================
# 4. EXPLORATORY VISUAL INSPECTION OF THE RESPONSE
# ==============================================================================
# A two-panel figure showing (a) the raw AOH distribution and (b) the
# log-transformed AOH distribution helps the reader verify that the chosen
# transformation effectively reshapes the response toward approximate normality.
message("\n[Diagnostics] Visual inspection of the response variable...")

plot_resp_raw <- ggplot(data_p3, aes(x = AOH_ha)) +
  geom_histogram(bins = 30, fill = "grey60", color = "black", alpha = 0.85) +
  labs(title = "(a) Raw OHE distribution",
       x = "OHE (ha)", y = "Frequency") +
  theme_classic(base_size = 10, base_family = "sans") +
  theme(plot.title = element_text(face = "bold", size = 11, hjust = 0))

plot_resp_log <- ggplot(data_p3, aes(x = log_AOH)) +
  geom_histogram(bins = 30, fill = "#009E73", color = "black", alpha = 0.85) +
  labs(title = "(b) log-transformed OHE distribution",
       x = "log(OHE)", y = "Frequency") +
  theme_classic(base_size = 10, base_family = "sans") +
  theme(plot.title = element_text(face = "bold", size = 11, hjust = 0))

plot_response_inspection <- plot_resp_raw + plot_resp_log + plot_layout(ncol = 2)
save_publication_plot(plot_response_inspection,
                      "Figure_PIII_01_Response_Transformation.tiff",
                      width_mm = 174, height_mm = 90)

# ==============================================================================
# 5. CHOOSING THE MODEL: SIX SPECIFICATIONS, THREE CRITERIA
# ==============================================================================
# chain of conditional rules: fit the OLS, test its assumptions, and fall back to
# another family if a test failed. That design has two defects. It never compares
# the families on equal terms, and it lets one marginal diagnostic settle a
# structural question. This section replaces the chain with an explicit
# comparison.
#
# THE ERROR THIS SECTION AVOIDS
# A model fitted to log(AOH) and a model fitted to AOH have likelihoods defined
# with respect to DIFFERENT response variables, so their information criteria are
# not comparable as they stand. The likelihood of a log-scale model is placed on
# the response scale by subtracting the Jacobian of the transformation,
# sum(log(y_i)). Without that correction the log-scale models appear about 1400
# AICc units better and the comparison is an artefact.
#
# THE THREE CRITERIA, in the order in which they are applied
#   1. Residual diagnostics. A specification whose residuals fail is not
#      considered, whatever its AICc. This is a gate, not a tie-breaker.
#   2. AICc weight on the common response scale, which measures relative support
#      given the data, and the number of parameters, which is the parsimony term
#      already inside the AICc.
#   3. Out-of-sample error in hectares, by exact leave-one-out cross validation,
#      which assumes no distribution and settles ties that the AICc leaves open.
# The nature of the response is the standing constraint behind all three: the
# occupied habitat envelope is continuous, strictly positive and right-skewed,
# and its dispersion grows with its mean.
# ==============================================================================
message("\n[Model choice] Comparing six specifications on a common response scale...")

options(na.action = "na.fail")   # required by MuMIn::dredge later on

ADD_RHS <- "PD_Forest_z + Mean_FRAC_Forest_z + Mean_ENN_Forest_m_z + PROX_MN_Forest_z + PLAND_Forest_z"
INT_RHS <- paste(ADD_RHS, "+ PROX_MN_Forest_z:PLAND_Forest_z")
JAC     <- sum(log(data_p3$AOH_ha))
n_obs   <- nrow(data_p3)
message(sprintf("   n = %d | Jacobian sum(log(AOH)) = %.2f", n_obs, JAC))

aicc_manual <- function(ll, k, n) -2 * ll + 2 * k + (2 * k * (k + 1)) / (n - k - 1)

# BIC is reported alongside AICc, not to replace it. Brewer, Butler and Cooksley
# (2016) show that the relative performance of AIC, AICc and BIC depends on the
# unobserved heterogeneity in the data, and these data carry heterogeneity that
# no predictor captures: 67 units from 9 taxa. Where the two criteria agree the
# choice is robust to that dependence; where they disagree the disagreement is
# itself the result, and is reported rather than resolved by preference.
bic_manual <- function(ll, k, n) -2 * ll + k * log(n)

ajustar <- function(familia, rhs) {
  if (familia == "OLS")   return(lm(as.formula(paste("log_AOH ~", rhs)), data = data_p3))
  if (familia == "LMM")   return(lme4::lmer(as.formula(paste("log_AOH ~", rhs, "+ (1 | GENUS)")),
                                            data = data_p3, REML = FALSE))
  if (familia == "Gamma") return(glm(as.formula(paste("AOH_ha ~", rhs)),
                                     data = data_p3, family = Gamma(link = "log")))
}

especificacoes <- expand.grid(
  familia  = c("OLS", "LMM", "Gamma"),
  estrutura = c("additive", "interaction"),
  stringsAsFactors = FALSE
)

modelos <- list(); linhas <- list()
for (i in seq_len(nrow(especificacoes))) {
  fam <- especificacoes$familia[i]; est <- especificacoes$estrutura[i]
  rhs <- if (est == "additive") ADD_RHS else INT_RHS
  m   <- ajustar(fam, rhs)
  rotulo <- sprintf("%s, %s", fam, est)
  modelos[[rotulo]] <- m

  # number of estimated parameters, including the variance components
  k <- switch(fam,
              OLS   = length(coef(m)) + 1,
              LMM   = length(lme4::fixef(m)) + 2,
              Gamma = length(coef(m)) + 1)
  # likelihood on the response scale
  ll <- if (fam == "Gamma") as.numeric(logLik(m)) else as.numeric(logLik(m)) - JAC

  # residual diagnostics by simulation, which apply to every family alike
  set.seed(GLOBAL_SEED)
  sim <- tryCatch(DHARMa::simulateResiduals(m, n = 1000, plot = FALSE),
                  error = function(e) NULL)
  if (is.null(sim)) {
    p_unif <- NA_real_; p_disp <- NA_real_; p_out <- NA_real_
  } else {
    p_unif <- DHARMa::testUniformity(sim, plot = FALSE)$p.value
    p_disp <- DHARMa::testDispersion(sim, plot = FALSE)$p.value
    p_out  <- DHARMa::testOutliers(sim, plot = FALSE)$p.value
  }
  linhas[[rotulo]] <- data.frame(
    Model = rotulo, Family = fam, Structure = est, k = k,
    logLik_response = round(ll, 3),
    AICc = round(aicc_manual(ll, k, n_obs), 3),
    BIC  = round(bic_manual(ll, k, n_obs), 3),
    DHARMa_uniformity = round(p_unif, 4),
    DHARMa_dispersion = round(p_disp, 4),
    DHARMa_outliers   = round(p_out, 4),
    stringsAsFactors = FALSE
  )
}
tab_comp <- do.call(rbind, linhas)
tab_comp$Delta_AICc <- round(tab_comp$AICc - min(tab_comp$AICc), 3)
tab_comp$Delta_BIC  <- round(tab_comp$BIC  - min(tab_comp$BIC),  3)
tab_comp$Weight     <- round(exp(-0.5 * tab_comp$Delta_AICc) /
                               sum(exp(-0.5 * tab_comp$Delta_AICc)), 4)
# a specification passes the gate when no simulation-based test rejects
tab_comp$Residuals_OK <- with(tab_comp,
  (is.na(DHARMa_uniformity) | DHARMa_uniformity > 0.05) &
  (is.na(DHARMa_dispersion) | DHARMa_dispersion > 0.05) &
  (is.na(DHARMa_outliers)   | DHARMa_outliers   > 0.05))

# ------------------------------------------------------------------------------
# 5.1 Out-of-sample error, in hectares
# ------------------------------------------------------------------------------
# Computed only for the families whose predictions are defined on the hectare
# scale without conditioning on a random effect, which is what a comparison of
# predictive accuracy across families requires.
# ------------------------------------------------------------------------------
# Leave-one-out is used rather than random k-fold. Yates, Aandahl, Richards and
# Brook (2022) recommend exact or approximate leave-one-out to minimize bias, or
# k-fold with a bias correction when k < 10; the earlier ten-fold scheme had
# neither. Leave-one-out is exact here because the models are cheap and n is 67.
#
# A second scheme, leave-one-genus-out, is reported beside it. The 67 units come
# from 9 taxa and units of the same taxon share species-level traits, so a random
# split places relatives in the training and the test set at once. Comparing the
# two quantifies how much that matters instead of assuming it away.
# ------------------------------------------------------------------------------
predict_ha <- function(familia, rhs, tr, te) {
  if (familia == "Gamma") {
    predict(glm(as.formula(paste("AOH_ha ~", rhs)), data = tr,
                family = Gamma(link = "log")), newdata = te, type = "response")
  } else {
    mm <- lm(as.formula(paste("log_AOH ~", rhs)), data = tr)
    # log-normal back-transformation, so that the two families are compared on
    # the hectare scale rather than on scales that are not the same quantity
    exp(predict(mm, newdata = te) + summary(mm)$sigma^2 / 2)
  }
}

# NOTE ON WHAT THESE TWO FUNCTIONS RETURN, added 22 Aug 2026.
# Both accumulate SQUARED errors on the hectare scale and return their square
# root. The quantity is therefore the ROOT MEAN SQUARED ERROR, and the exported
# columns are named CV_RMSE_LOO_ha and CV_RMSE_LGO_ha accordingly. An earlier
# draft of the manuscript described these values as a mean absolute error, which
# they are not; the two differ whenever the error distribution is skewed, and
# here it is dominated by one unit of 387,619 ha. Report them as RMSE.
cv_rmse_loo <- function(familia, rhs) {
  err <- numeric(0)
  for (i in seq_len(n_obs)) {
    tr <- data_p3[-i, ]; te <- data_p3[i, , drop = FALSE]
    p  <- tryCatch(predict_ha(familia, rhs, tr, te), error = function(e) NA_real_)
    err <- c(err, (te$AOH_ha - p)^2)
  }
  sqrt(mean(err, na.rm = TRUE))
}

cv_rmse_lgo <- function(familia, rhs) {
  err <- numeric(0)
  for (g in unique(data_p3$GENUS)) {
    tr <- data_p3[data_p3$GENUS != g, ]; te <- data_p3[data_p3$GENUS == g, ]
    p  <- tryCatch(predict_ha(familia, rhs, tr, te),
                   error = function(e) rep(NA_real_, nrow(te)))
    err <- c(err, (te$AOH_ha - p)^2)
  }
  sqrt(mean(err, na.rm = TRUE))
}

tab_comp$CV_RMSE_LOO_ha <- NA_real_
tab_comp$CV_RMSE_LGO_ha <- NA_real_
for (rotulo in tab_comp$Model) {
  fam <- tab_comp$Family[tab_comp$Model == rotulo]
  est <- tab_comp$Structure[tab_comp$Model == rotulo]
  if (fam == "LMM") next
  rhs <- if (est == "additive") ADD_RHS else INT_RHS
  tab_comp$CV_RMSE_LOO_ha[tab_comp$Model == rotulo] <- round(cv_rmse_loo(fam, rhs), 1)
  tab_comp$CV_RMSE_LGO_ha[tab_comp$Model == rotulo] <- round(cv_rmse_lgo(fam, rhs), 1)
}
# kept under the old name so that downstream code and the tie-break rule below
# continue to read the primary criterion
tab_comp$CV_RMSE_ha <- tab_comp$CV_RMSE_LOO_ha

tab_comp <- tab_comp[order(tab_comp$AICc), ]
write.csv(tab_comp, file.path(output_path, "Table_S_PartIII_ModelComparison.csv"),
          row.names = FALSE)
message("\n   Comparison of the six specifications:")
print(tab_comp[, c("Model", "k", "AICc", "Delta_AICc", "Weight", "Delta_BIC",
                   "Residuals_OK", "CV_RMSE_LOO_ha", "CV_RMSE_LGO_ha")], row.names = FALSE)

if (tab_comp$Model[which.min(tab_comp$AICc)] != tab_comp$Model[which.min(tab_comp$BIC)]) {
  message("   NOTE. AICc and BIC point to different specifications:")
  message(sprintf("         AICc favours %s; BIC favours %s.",
                  tab_comp$Model[which.min(tab_comp$AICc)],
                  tab_comp$Model[which.min(tab_comp$BIC)]))
  message("         Report the disagreement; it reflects the unobserved heterogeneity")
  message("         among taxa (Brewer et al. 2016) and is not resolved by preference.")
} else {
  message("   AICc and BIC agree on the same specification.")
}

# ------------------------------------------------------------------------------
# 5.2 The decision
# ------------------------------------------------------------------------------
# Applied in the stated order. Among the specifications whose residuals pass,
# those within two AICc units of the best are treated as indistinguishable by
# that criterion, following Burnham & Anderson (2002); the tie is then settled by
# the out-of-sample error, and, where that is also close, by parsimony.
elegiveis <- tab_comp[tab_comp$Residuals_OK, , drop = FALSE]
if (nrow(elegiveis) == 0) {
  stop("No specification passed the residual diagnostics. Inspect Table_S_PartIII_ModelComparison.csv.")
}
empatados <- elegiveis[elegiveis$Delta_AICc <= 2, , drop = FALSE]
if (all(is.na(empatados$CV_RMSE_ha))) {
  escolhido <- empatados$Model[which.min(empatados$k)]
  criterio  <- "fewest parameters among the models tied on AICc"
} else {
  escolhido <- empatados$Model[which.min(empatados$CV_RMSE_ha)]
  criterio  <- "lowest cross-validated error among the models tied on AICc"
}
global_model       <- modelos[[escolhido]]
model_family_label <- escolhido

message(sprintf("\n[Decision] Retained specification: %s", escolhido))
message(sprintf("            Criterion: %s.", criterio))
message(sprintf("            Delta AICc = %.2f, weight = %.3f, k = %d.",
                tab_comp$Delta_AICc[tab_comp$Model == escolhido],
                tab_comp$Weight[tab_comp$Model == escolhido],
                tab_comp$k[tab_comp$Model == escolhido]))
message("            The structure of the response is the standing constraint:")
message("            AOH is continuous, strictly positive and right-skewed, and")
message("            its dispersion grows with its mean.")

formula_global <- formula(global_model)

# Collinearity within the retained predictor pool. VIF is computed on the linear
# predictor, which is the quantity the diagnostic refers to in any of the three
# families. Values above 5 for the terms that make up the product are expected
# and are not an artefact of the product itself, since the predictors were
# standardized before it was formed (Aiken & West 1991; Schielzeth 2010).
vif_partIII <- car::vif(lm(as.formula(paste("log_AOH ~", INT_RHS)), data = data_p3))
write.csv(data.frame(Variable = names(vif_partIII), VIF = round(vif_partIII, 3)),
          file.path(output_path, "Table_S_PartIII_VIF.csv"), row.names = FALSE)
message(sprintf("   Maximum VIF among the terms of the retained model: %.2f", max(vif_partIII)))

# Diagnostics of the retained model, exported on their own.
diag_final <- run_residual_diagnostics(global_model, label = model_family_label)
write.csv(diag_final,
          file.path(output_path, "Table_S_PartIII_Diagnostics.csv"), row.names = FALSE)

# ------------------------------------------------------------------------------
# 5.3 Taxonomic dependence, reported rather than decided by
# ------------------------------------------------------------------------------
# The 67 units come from 9 taxa, so units of the same taxon are not independent.
# The mixed specification is in the comparison above and is not favoured by it;
# the intraclass correlation is reported because a reader will ask for it, not
# because a threshold on it decides anything.
lmm_ref <- modelos[["LMM, interaction"]]
vc      <- as.data.frame(lme4::VarCorr(lmm_ref))
icc_su  <- vc$vcov[vc$grp == "GENUS"] / sum(vc$vcov)
write.csv(
  data.frame(Test  = c("ICC (genus), from the mixed specification",
                       "Delta AICc, retained model minus mixed model",
                       "Weight of the mixed specification"),
             Value = round(c(icc_su,
                             tab_comp$AICc[tab_comp$Model == escolhido] -
                               tab_comp$AICc[tab_comp$Model == "LMM, interaction"],
                             tab_comp$Weight[tab_comp$Model == "LMM, interaction"]), 4)),
  file.path(output_path, "Table_S_PartIII_DependenceDiagnostic.csv"), row.names = FALSE)
message(sprintf("   ICC of the genus in the mixed specification: %.3f", icc_su))

# Global-model summary (single-row record).
glob_summary <- data.frame(
  Family            = model_family_label,
  N_obs             = n_obs,
  N_predictors      = length(all.vars(formula_global)) - 1,
  R2_or_DevExpl     = round(explained_variance(global_model), 4),
  AICc              = round(tab_comp$AICc[tab_comp$Model == escolhido], 3),
  stringsAsFactors  = FALSE
)
write.csv(glob_summary,
          file.path(output_path, "Table_S_PartIII_GlobalModelSummary.csv"),
          row.names = FALSE)
message(sprintf("   Explained variance / deviance of the retained model: %.3f",
                glob_summary$R2_or_DevExpl))


# ==============================================================================
# 7.5 THE CONFIRMATORY MODEL: ONE SPECIFICATION, FITTED ONCE
# ==============================================================================
# This is the inferential result of the script, and it is deliberately separated
# from everything that follows.
#
# Yates et al. (2023) state the problem plainly: inference on parameter estimates
# is biased when it is preceded by model selection, and valid inference requires
# either a carefully specified single model or technical post-selection
# adjustments. Tredennick et al. (2021) make the same distinction by goal: the
# procedure that suits prediction does not suit inference, and the confusion in
# the literature comes from not declaring the goal first.
#
# The goal of this study is inference. The candidate set was defined a priori,
# and the full a priori model is therefore fitted once, without selection, and
# its coefficients are the estimates that carry a confidence interval and a p
# value. Section 8 explores the subsets, and what it produces is a description of
# selection uncertainty, not a second set of estimates to be tested.
#
# Nothing here is chosen after seeing the data. The family was chosen in Section
# 5 by residual behavior and out-of-sample error, both of which are properties
# of the specification rather than of any coefficient.
# ==============================================================================
message("\n[Confirmatory] Fitting the full a priori model, without selection...")

confirm_model <- global_model
cf_conf <- summary(confirm_model)$coefficients
ci_conf <- confint.default(confirm_model)

confirm_table <- data.frame(
  Predictor = rownames(cf_conf),
  Estimate  = round(cf_conf[, 1], 4),
  SE        = round(cf_conf[, 2], 4),
  CI_lower  = round(ci_conf[, 1], 4),
  CI_upper  = round(ci_conf[, 2], 4),
  Statistic = round(cf_conf[, 3], 4),
  P_value   = signif(cf_conf[, 4], 4),
  Source    = "Full a priori model, no selection",
  row.names = NULL, stringsAsFactors = FALSE
)
write.csv(confirm_table,
          file.path(output_path, "Table_S_PartIII_ConfirmatoryModel.csv"),
          row.names = FALSE)
message("   Coefficients of the confirmatory model:")
print(confirm_table[, c("Predictor", "Estimate", "SE", "CI_lower", "CI_upper", "P_value")],
      row.names = FALSE)


# ==============================================================================
# 7.6 SPATIAL SENSITIVITY OF THE CONFIRMATORY MODEL
# ==============================================================================
# Diniz-Filho, Rangel and Bini (2008) show that an information criterion computed
# on a non-spatial regression is sensitive to spatial autocorrelation and yields
# unstable, overfitted minimum adequate models. Script 02 detected spatial
# structure in these data: the global Moran eigenvector model is significant and
# two eigenvectors were retained.
#
# The selection above is non-spatial. Rather than leave that as an unquantified
# caveat, the confirmatory model is refitted with the retained eigenvectors as
# covariates and the coefficients are compared. If they hold, the non-spatial
# selection did not distort them; if they move, the manuscript must say so.
# ==============================================================================
message("\n[Spatial] Refitting the confirmatory model with the retained Moran eigenvectors...")

mem_file <- file.path(output_path, "PartII", "TableII_2b_MEM_Vectors.csv")

if (file.exists(mem_file)) {
  mem_df <- read_csv(mem_file, show_col_types = FALSE) %>%
    mutate(UA_ID = as.character(UA_ID), SPECIES = as.character(SPECIES))
  mem_cols <- setdiff(names(mem_df), c("SPECIES", "UA_ID"))

  d_sp <- data_p3 %>%
    mutate(UA_ID = as.character(UA_ID), SPECIES = as.character(SPECIES)) %>%
    left_join(mem_df, by = c("SPECIES", "UA_ID"))

  if (anyNA(d_sp[, mem_cols])) {
    message("   Eigenvectors did not join to every unit; the check is skipped.")
  } else {
    rhs_sp  <- paste(INT_RHS, "+", paste(mem_cols, collapse = " + "))
    m_sp    <- glm(as.formula(paste("AOH_ha ~", rhs_sp)), data = d_sp,
                   family = Gamma(link = "log"))
    cf_sp   <- summary(m_sp)$coefficients
    shared  <- intersect(rownames(cf_conf), rownames(cf_sp))

    spatial_table <- data.frame(
      Predictor          = shared,
      Beta_non_spatial   = round(cf_conf[shared, 1], 4),
      Beta_with_MEMs     = round(cf_sp[shared, 1], 4),
      Absolute_change    = round(abs(cf_sp[shared, 1] - cf_conf[shared, 1]), 4),
      P_non_spatial      = signif(cf_conf[shared, 4], 4),
      P_with_MEMs        = signif(cf_sp[shared, 4], 4),
      stringsAsFactors   = FALSE, row.names = NULL
    )
    write.csv(spatial_table,
              file.path(output_path, "Table_S_PartIII_SpatialSensitivity.csv"),
              row.names = FALSE)
    print(spatial_table, row.names = FALSE)

    flips <- with(spatial_table,
                  sign(Beta_non_spatial) != sign(Beta_with_MEMs) &
                    Predictor != "(Intercept)")
    if (any(flips)) {
      message("   WARNING. At least one coefficient changes sign when the spatial")
      message("   structure is included. Report the spatial specification.")
    } else {
      message("   No coefficient changes sign. The non-spatial selection did not")
      message("   distort the estimates, and this is the sentence the Methods needs.")
    }
  }
} else {
  message("   TableII_2b_MEM_Vectors.csv not found. Run script 02 first; without it")
  message("   the spatial sensitivity of Diniz-Filho et al. (2008) cannot be checked.")
}


# ==============================================================================
# 8. EXPLORATORY MULTIMODEL ANALYSIS (AICc, dredge, model averaging)
# ==============================================================================
# All 2^5 = 32 candidate subsets that respect marginal are ranked by AICc
# with MuMIn::dredge: five predictors plus the a priori product, the product
# admitted only where both constituent terms are present.
# The competitive set is defined by Delta AICc < 2 (Burnham & Anderson 2002).
#
# When two or more models share Delta AICc < 2, model averaging is required:
# ignoring the second-best model would discard genuine selection uncertainty,
# while full averaging would bias coefficients of predictors that are absent
# from some competitive models toward zero. The conventional solution is the
# subset (conditional) average (Burnham & Anderson 2002, Grueber et al. 2011):
# each parameter is averaged only across the models in which it actually
# occurs. This isolates the "pure" effect of each predictor while preserving
# the selection uncertainty in the variable-importance metric (sum of Akaike
# weights).
#
# In this dataset, two competitive models (Delta AICc < 2) are typically
# retained: the difference between them is the inclusion or exclusion of
# Mean_ENN_Forest_m. The subset average produces standardized coefficients
# that reflect the average effect of each predictor where it matters, and
# acknowledges the model-selection uncertainty for forest isolation (ENN_MN).
# ==============================================================================
# ------------------------------------------------------------------------------
# WHAT THIS SECTION IS, AND WHAT IT IS NOT
#
# It is a description of selection uncertainty across the a priori candidate set.
# It is NOT the inferential result: the coefficients below were produced after a
# search over 16 subsets, so their standard errors and p values are conditional
# on a selection step and are biased as inference (Yates et al. 2023). The
# estimates that the manuscript interprets are those of Section 7.5.
#
# Both averaging schemes are reported, because they answer different questions.
# The subset, or conditional, average takes each parameter over the models that
# contain it, and describes the effect where it appears. The full average enters
# a zero for the models that omit the parameter, and is the more conservative
# summary when the question is how much support the predictor has across the
# whole set. Reporting only one of them hides the difference between "strong
# where it appears" and "appears often".
# ------------------------------------------------------------------------------
message("\n[Exploratory] AICc-based multimodel analysis of the candidate set...")

set.seed(GLOBAL_SEED)
# subset = dc(...) enforces marginality: the product is only ever fitted in
# models that already contain both of its constituent terms.
# dredge runs on the retained specification, so the candidate set is explored
# within the family already chosen in Section 5 rather than across families.
# subset = dc(...) enforces marginality: the product is only fitted in models
# that already contain both of its constituent terms.
dredge_table <- MuMIn::dredge(
  global_model, rank = "AICc", trace = FALSE,
  subset = dc(PROX_MN_Forest_z, PLAND_Forest_z, `PROX_MN_Forest_z:PLAND_Forest_z`)
)

# Export the full dredge table (all subsets that respect marginality).
dredge_df <- as.data.frame(dredge_table) %>%
  mutate(across(where(is.numeric), ~round(.x, 4)))
write.csv(dredge_df,
          file.path(output_path, "Table_S_PartIII_Dredge_FullTable.csv"),
          row.names = FALSE)

# Identify competitive set (Delta AICc < 2; Burnham & Anderson 2002).
competitive <- MuMIn::get.models(dredge_table, subset = delta < 2)
n_competitive <- length(competitive)
message(sprintf("   Competitive set (Delta AICc < 2): %d model(s).",
                n_competitive))

# Multimodel inference logic. Subset (conditional) average is used by
# extract_coefs() when n_competitive > 1; see Section 8 header for the
# rationale.
if (n_competitive == 1) {
  final_model <- competitive[[1]]
  averaged    <- FALSE
  inf_label   <- "Single most-parsimonious model retained"
} else {
  final_model <- MuMIn::model.avg(competitive, fit = TRUE)
  averaged    <- TRUE
  inf_label   <- sprintf(
    "Conditional (subset) average across %d competitive models (Delta AICc < 2)",
    n_competitive
  )
}
message(sprintf("   Inference: %s.", inf_label))

# Variable importance via the sum of Akaike weights across the full dredge.
sw_vec <- MuMIn::sw(dredge_table)
df_importance <- data.frame(
  Predictor   = names(sw_vec),
  Importance  = round(as.numeric(sw_vec), 4),
  stringsAsFactors = FALSE
) %>%
  arrange(desc(Importance))
write.csv(df_importance,
          file.path(output_path, "Table_S_PartIII_VariableImportance.csv"),
          row.names = FALSE)

# Coefficient table for the final model (single or averaged).
coef_df <- extract_coefs(final_model, label = inf_label) %>%
  mutate(across(where(is.numeric), ~round(.x, 4)))

# Full (zero-method) average alongside the conditional one, with the revised
# variance estimator of Burnham and Anderson (2004), which MuMIn implements as
# revised.var and which widens the interval to acknowledge selection uncertainty.
if (averaged) {
  avg_rev <- MuMIn::model.avg(competitive, revised.var = TRUE)
  ct_full <- MuMIn::coefTable(avg_rev, full = TRUE)
  ct_cond <- MuMIn::coefTable(avg_rev, full = FALSE)
  ci_full <- confint(avg_rev, full = TRUE)
  sw_all  <- MuMIn::sw(dredge_table)

  both <- data.frame(
    Predictor      = rownames(ct_full),
    Beta_full      = round(ct_full[, 1], 4),
    SE_full        = round(ct_full[, 2], 4),
    CI_lower_full  = round(ci_full[, 1], 4),
    CI_upper_full  = round(ci_full[, 2], 4),
    Beta_conditional = round(ct_cond[match(rownames(ct_full), rownames(ct_cond)), 1], 4),
    SE_conditional   = round(ct_cond[match(rownames(ct_full), rownames(ct_cond)), 2], 4),
    stringsAsFactors = FALSE, row.names = NULL
  )
  both$Sum_of_weights <- round(sw_all[match(both$Predictor, names(sw_all))], 4)
  write.csv(both,
            file.path(output_path, "Table_S_PartIII_Averaging_FullVsConditional.csv"),
            row.names = FALSE)
  message("\n   Full and conditional averages, with the sum of Akaike weights:")
  print(both, row.names = FALSE)
  message("   These describe selection uncertainty. They are not the inferential")
  message("   estimates; those are in Table_S_PartIII_ConfirmatoryModel.csv.")
}
write.csv(coef_df,
          file.path(output_path, "Table_S_PartIII_FinalCoefficients.csv"),
          row.names = FALSE)

# Compact comparison table of the competitive models for reporting in the main
# results.
comp_df <- as.data.frame(subset(dredge_table, delta < 2))
metric_cols <- c("df", "logLik", "AICc", "delta", "weight")
var_cols    <- setdiff(names(comp_df), metric_cols)
formulas <- vapply(seq_len(nrow(comp_df)), function(i) {
  vp <- var_cols[!is.na(comp_df[i, var_cols])]
  vp <- vp[vp != "(Intercept)"]
  if (length(vp) == 0) "~ Null Model ~" else paste(vp, collapse = " + ")
}, character(1))
comp_clean <- data.frame(
  Rank               = seq_len(nrow(comp_df)),
  Selected_Variables = formulas,
  df                 = comp_df$df,
  LogLik             = round(comp_df$logLik, 3),
  AICc               = round(comp_df$AICc,   3),
  Delta_AICc         = round(comp_df$delta,  3),
  Akaike_Weight      = round(comp_df$weight, 3)
)
write.csv(comp_clean,
          file.path(output_path, "Table_S_PartIII_CompetitiveModels.csv"),
          row.names = FALSE)

# ==============================================================================
# 8.5 COLLINEARITY STABILITY CHECK: PD_Forest x Mean_FRAC_Forest
# ==============================================================================
# PD_Forest and Mean_FRAC_Forest are moderately correlated and enter the retained
# model with OPPOSITE signs (PD negative, FRAC positive). The pairwise Pearson
# correlation is computed below on the 67 Part III records and written to
# Table_S_PartIII_*.csv. On the earlier 65-record dataset it was about 0.66,
# below the 0.7 screening threshold of Dormann et al. (2013), with multivariable
# VIF well below 5 (max 2.2). Quote the exported value, not this comment.
# Collinearity is moderate rather than severe as long as the printed r stays
# below 0.7 and the printed VIF below 5; check both after each run.
# The opposite signs are still worth documenting, because two correlated
# predictors with opposite-sign effects can reinforce each other in the joint
# model (mutual suppression). To check this, two reduced OLS models are fitted,
# one dropping PD_Forest and one dropping Mean_FRAC_Forest, and the coefficients
# are compared with the full model. If a surviving predictor keeps its sign and a
# comparable magnitude, the joint interpretation is defensible; if a coefficient
# collapses, the Discussion must read the pair with explicit caution.
#
# Why Pearson here. Collinearity is reported with Pearson because the model that
# the check protects (the OLS, and the VIF) is linear, and VIF measures the
# LINEAR dependence among predictors. The pairwise Pearson correlation is the
# matching diagnostic for that linear dependence (Dormann et al. 2013), and it
# does not require the predictors to be normally distributed; normality concerns
# only the parametric significance test of r, not r as a measure of linear
# covariation. Even so, a Pearson value can be distorted by skew or outliers, so
# the Spearman (rank) correlation is reported alongside it as a robustness
# cross-check (Zuur et al. 2010). When the two agree, the linear reading is safe.
# Spearman is the primary statistic in the descriptive per-genus heatmaps
# (Sections 10.5 and 10.6 companions), where samples are small and associations
# need not be linear; there, a distribution-free, outlier-robust measure is the
# better default.
# This block reports a Supplementary table; it does not change the main model.
# ==============================================================================
message("\n[Stability] Collinearity stability check (PD vs FRAC)...")

# Pairwise correlation, reported for the record. Pearson is the collinearity
# diagnostic that matches the linear model and VIF; Spearman is the
# distribution-free cross-check (see the section header).
r_pd_frac          <- cor(data_p3$PD_Forest, data_p3$Mean_FRAC_Forest,
                          method = "pearson", use = "complete.obs")
rho_pd_frac        <- cor(data_p3$PD_Forest, data_p3$Mean_FRAC_Forest,
                          method = "spearman", use = "complete.obs")
message(sprintf("   r(PD_Forest, Mean_FRAC_Forest): Pearson = %.3f, Spearman = %.3f",
                r_pd_frac, rho_pd_frac))

# Full four-predictor OLS (reference for the comparison).
ajusta_est <- function(rhs) {
  if (inherits(global_model, "glm")) {
    glm(as.formula(paste("AOH_ha ~", rhs)), data = data_p3, family = Gamma(link = "log"))
  } else {
    lm(as.formula(paste("log_AOH ~", rhs)), data = data_p3)
  }
}
# The stability check is run in the retained family so that the coefficients it
# reports are on the same scale as Table 3.
lm_full_4    <- ajusta_est(INT_RHS)
lm_drop_pd   <- ajusta_est(sub("PD_Forest_z \\+ ", "", INT_RHS))
lm_drop_frac <- ajusta_est(sub("Mean_FRAC_Forest_z \\+ ", "", INT_RHS))

# Helper: extract one predictor's coefficient from a model, NA if absent.
# `term` is the first argument so it can be iterated with vapply() while the
# model is passed by name.
get_beta <- function(term, model) {
  cf <- coef(model)
  if (term %in% names(cf)) as.numeric(cf[term]) else NA_real_
}

stability_terms <- c("PD_Forest_z", "Mean_FRAC_Forest_z",
                     "Mean_ENN_Forest_m_z", "PROX_MN_Forest_z")

stability_df <- data.frame(
  Predictor      = stability_terms,
  Beta_full      = round(vapply(stability_terms, get_beta,
                                numeric(1), model = lm_full_4), 4),
  Beta_drop_PD   = round(vapply(stability_terms, get_beta,
                                numeric(1), model = lm_drop_pd), 4),
  Beta_drop_FRAC = round(vapply(stability_terms, get_beta,
                                numeric(1), model = lm_drop_frac), 4),
  row.names = NULL
)
# Append the pairwise correlations as footnote rows for traceability: Pearson
# (the collinearity diagnostic that matches the linear model) and Spearman (the
# distribution-free cross-check). Close agreement supports the linear reading.
stability_df <- rbind(
  stability_df,
  data.frame(Predictor = "Pearson r(PD_Forest, Mean_FRAC_Forest)",
             Beta_full = round(r_pd_frac, 4),
             Beta_drop_PD = NA, Beta_drop_FRAC = NA),
  data.frame(Predictor = "Spearman rho(PD_Forest, Mean_FRAC_Forest)",
             Beta_full = round(rho_pd_frac, 4),
             Beta_drop_PD = NA, Beta_drop_FRAC = NA)
)
write.csv(stability_df,
          file.path(output_path, "Table_S_PartIII_StabilityCheck.csv"),
          row.names = FALSE)
message("   Stability check exported (full vs PD-dropped vs FRAC-dropped).")

# ==============================================================================
# 8.6 SENSITIVITY MODEL: ADD MATRIX PROXIMITY (PROX_MN_Agropecuaria)
# ==============================================================================
# The farming-matrix proximity term is the strongest term in the Part I
# RDA but was deliberately excluded from the main Part III model so that the
# response (AOH inside the remnant) is explained by forest configuration only.
# Reviewers commonly ask for the matrix effect on AOH itself. This block fits
# the five-predictor OLS that adds PROX_MN_Agropecuaria, runs the same dredge
# selection, and exports the coefficients of the best model as a sensitivity
# analysis for the Supplementary Material. The block self-skips, with a clear
# message, if the matrix predictor is not present in the raw data.
# ==============================================================================
message("\n[Sensitivity] Matrix-proximity sensitivity model...")

matrix_pred <- "PROX_MN_Agropecuaria"
if (matrix_pred %in% names(data_raw)) {
  
  # Before v7 this drop_na() removed A_guariba UA 6 and B_arachnoides UA 1,
  # whose PROX_MN_Agropecuaria cell was empty, so the sensitivity model was
  # fitted on 65 observations and compared against a retained model fitted on 67.
  # A sensitivity check is only informative when the two models share a sample.
  data_sens <- data_raw %>%
    select(UA_ID, GENUS, AOH_ha, all_of(predictors_part3), all_of(matrix_pred)) %>%
    drop_na(all_of(c(predictors_part3, matrix_pred))) %>%
    mutate(AOH_ha  = ifelse(AOH_ha == 0, 0.001, AOH_ha),
           log_AOH = log(AOH_ha)) %>%
    mutate(across(all_of(c(predictors_part3, matrix_pred)),
                  ~as.numeric(scale(.x)), .names = "{.col}_z"))
  
  # The sensitivity model is fitted in the SAME family as the retained model,
  # otherwise its coefficients would not be comparable with Table 3.
  if (nrow(data_sens) != n_obs) {
    stop(sprintf(paste0("The sensitivity model would use %d observations while ",
                        "the retained model uses %d. The two must share a ",
                        "sample to be comparable. Check PROX_MN_Agropecuaria in ",
                        "Data_Raw_WithCoords.csv."),
                 nrow(data_sens), n_obs))
  }
  message(sprintf("   Sensitivity model fitted on %d observations, as the retained model.",
                  nrow(data_sens)))

  rhs_sens <- paste(INT_RHS, "+ PROX_MN_Agropecuaria_z")
  lm_sens <- if (inherits(global_model, "glm")) {
    glm(as.formula(paste("AOH_ha ~", rhs_sens)), data = data_sens,
        family = Gamma(link = "log"))
  } else {
    lm(as.formula(paste("log_AOH ~", rhs_sens)), data = data_sens)
  }
  
  # VIF on the five-predictor pool (the matrix term may reintroduce collinearity).
  vif_sens <- car::vif(lm_sens)
  message(sprintf("   Max VIF in 5-predictor sensitivity pool: %.2f",
                  max(vif_sens)))
  
  # Same AICc dredge on the larger model.
  set.seed(GLOBAL_SEED)
  dredge_sens <- MuMIn::dredge(lm_sens, rank = "AICc", trace = FALSE)
  best_sens   <- MuMIn::get.models(dredge_sens, subset = 1)[[1]]
  
  coef_sens <- extract_coefs(best_sens,
                             label = "Sensitivity: best 5-predictor subset") %>%
    mutate(across(where(is.numeric), ~round(.x, 4)))
  # Attach the matrix-term VIF as a trailing diagnostic column where available.
  coef_sens$VIF_full_model <- round(
    vif_sens[match(stringr::str_remove(coef_sens$Predictor, "_z$") %>%
                     paste0("_z"), names(vif_sens))], 3)
  write.csv(coef_sens,
            file.path(output_path, "Table_S_PartIII_Sensitivity_MatrixModel.csv"),
            row.names = FALSE)
  message("   Matrix-proximity sensitivity model exported.")
  
} else {
  message(sprintf("   '%s' not found in the input data; sensitivity model skipped.",
                  matrix_pred))
}

# ==============================================================================
# 8.7 DESCRIPTIVE CHECK: FRAC vs FOREST PATCH AREA
# ==============================================================================
# This block checks whether shape complexity (Mean_FRAC_Forest) tracks forest
# mean patch area (Mean_AREA_Forest_ha), to inform the reading of the positive
# FRAC effect on AOH. In this dataset the two are NEGATIVELY correlated: on the
# earlier 65-record version, Pearson about -0.67 and Spearman about -0.63. The
# coefficients are recomputed below on the 67 records and exported to
# Table_S_PartIII_FRAC_Area_Correlation.csv; quote those. More complex-shaped
# patches tend to be SMALLER, not larger. The positive FRAC effect on AOH therefore cannot be
# attributed to patch extent. It must be read as an effect of shape complexity
# itself, for example more convoluted forest boundaries that may track rugged
# terrain where forest persists, and the Discussion should frame the mechanism as
# a hypothesis rather than a settled size effect. Both Pearson and Spearman are
# reported, because the marginal distributions can be skewed; close agreement
# between the two indicates the result is not an artefact of skew or outliers.
# Descriptive only; this block does not change the main model.
# ==============================================================================
message("\n[Descriptive] FRAC vs forest patch area correlation...")

area_pred <- "Mean_AREA_Forest_ha"
if (area_pred %in% names(data_raw)) {
  area_check <- data_raw %>%
    select(Mean_FRAC_Forest, all_of(area_pred)) %>%
    drop_na()
  r_frac_area_pearson  <- cor(area_check$Mean_FRAC_Forest,
                              area_check[[area_pred]], method = "pearson")
  r_frac_area_spearman <- cor(area_check$Mean_FRAC_Forest,
                              area_check[[area_pred]], method = "spearman")
  frac_area_df <- data.frame(
    Metric = c("Pearson r (FRAC vs patch area)",
               "Spearman rho (FRAC vs patch area)",
               "N"),
    Value  = round(c(r_frac_area_pearson, r_frac_area_spearman,
                     nrow(area_check)), 4)
  )
  write.csv(frac_area_df,
            file.path(output_path, "Table_S_PartIII_FRAC_Area_Correlation.csv"),
            row.names = FALSE)
  message(sprintf("   r(FRAC, patch area): Pearson = %.3f, Spearman = %.3f (n = %d).",
                  r_frac_area_pearson, r_frac_area_spearman, nrow(area_check)))
} else {
  message(sprintf("   '%s' not in input; FRAC-area check skipped.", area_pred))
}

# ==============================================================================
# 9. FIGURE: COEFFICIENT FOREST PLOT (standardized effects with 95% CI)
# ==============================================================================
# Figure S10 of the supplement. The panel is split in two, and the split is the
# point of the figure: three of the six terms of the a priori set are normalized
# by the response, or built on a term that is, and the permutation null of
# Part V reproduces all three (Table S28b). They are drawn hollow and in grey and
# carry no significance mark, because a star beside a coefficient the null
# reproduces invites the reading this paper rejects. The three terms free of the
# response keep the mark.
#
# REVISED 24 Aug 2026. Until then the figure sorted all six terms together and
# starred them all, which showed the a priori model as if every coefficient could
# be read.
# ==============================================================================
message("\n[Visualization] Coefficient forest plot...")

# The terms the response normalises. PD is NP divided by the unit area, PLAND is
# the class area divided by the unit area, and the product inherits PLAND.
coupled_terms <- c("PD_Forest_z", "PLAND_Forest_z",
                   "PLAND_Forest_z:PROX_MN_Forest_z",
                   "PROX_MN_Forest_z:PLAND_Forest_z")

GRP_FREE    <- "Free of the response: interpretable"
GRP_COUPLED <- "Normalised by the response: not interpretable"

forest_df <- coef_df %>%
  filter(Predictor != "(Intercept)") %>%
  mutate(
    Predictor_raw = strip_z(Predictor),
    Label         = label_predictors(Predictor_raw),
    Group         = ifelse(Predictor %in% coupled_terms, GRP_COUPLED, GRP_FREE),
    Significance  = ifelse(Predictor %in% coupled_terms, "",
                           case_when(P_value < 0.001 ~ "***",
                                     P_value < 0.01  ~ "**",
                                     P_value < 0.05  ~ "*",
                                     TRUE            ~ ""))
  ) %>%
  arrange(Group, Estimate) %>%
  mutate(Label = factor(Label, levels = Label),
         Group = factor(Group, levels = c(GRP_FREE, GRP_COUPLED)))

plot_forest <- ggplot(forest_df,
                      aes(x = Estimate, y = Label,
                          colour = Group, fill = Group)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
  geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper),
                 height = 0.18, linewidth = 0.7) +
  geom_point(size = 3.4, shape = 21, stroke = 0.9) +
  geom_text(aes(label = Significance,
                x = ifelse(Estimate >= 0, CI_upper + 0.04, CI_lower - 0.04)),
            colour = "black", size = 4.2, fontface = "bold",
            hjust = ifelse(forest_df$Estimate >= 0, 0, 1),
            show.legend = FALSE) +
  facet_wrap(~ Group, ncol = 1, scales = "free_y", strip.position = "top") +
  scale_colour_manual(values = setNames(c("#0072B2", "#8A8A8A"),
                                        c(GRP_FREE, GRP_COUPLED)), guide = "none") +
  scale_fill_manual(values = setNames(c("#0072B2", "white"),
                                      c(GRP_FREE, GRP_COUPLED)), guide = "none") +
  labs(x = expression("Standardised coefficient (" * beta * ", 95% CI)"),
       y = NULL) +
  theme_classic(base_size = 11, base_family = "sans") +
  theme(axis.text.y        = element_text(size = 9.5, face = "bold"),
        axis.title.x       = element_text(face = "bold", size = 10),
        strip.text         = element_text(size = 9, face = "bold", hjust = 0,
                                          margin = margin(3, 3, 3, 3)),
        strip.background   = element_rect(fill = "grey95", colour = NA),
        panel.grid.major.x = element_line(colour = "grey92", linewidth = 0.3),
        panel.spacing      = unit(6, "pt"))

save_publication_plot(plot_forest,
                      "Figure_PIII_02_Coefficient_ForestPlot.tiff",
                      width_mm = 174, height_mm = 108)

# ==============================================================================
# 10. FIGURE: VARIABLE IMPORTANCE (Akaike weights)
# ==============================================================================
message("\n[Visualization] Variable importance bar chart...")

imp_plot_df <- df_importance %>%
  mutate(
    Predictor_clean = strip_z(Predictor),
    Label = label_predictors(Predictor_clean)
  ) %>%
  arrange(Importance) %>%
  mutate(Label = factor(Label, levels = Label))

plot_importance <- ggplot(imp_plot_df,
                          aes(x = Importance, y = Label, fill = Importance)) +
  geom_col(color = "black", linewidth = 0.4, width = 0.7) +
  geom_text(aes(label = sprintf("%.2f", Importance)),
            hjust = -0.15, size = 3.4, fontface = "bold") +
  scale_fill_viridis_c(option = "mako", direction = -1, limits = c(0, 1),
                       name = "Akaike\nweight") +
  scale_x_continuous(limits = c(0, 1.1), breaks = seq(0, 1, 0.25)) +
  labs(x = "Sum of Akaike weights",
       y = NULL) +
  theme_classic(base_size = 11, base_family = "sans") +
  theme(plot.title = element_text(face = "bold", size = 11, hjust = 0.5),
        axis.text.y = element_text(size = 9.5, face = "bold"),
        legend.position = "right")

save_publication_plot(plot_importance,
                      "Figure_PIII_03_Variable_Importance.tiff",
                      width_mm = 174, height_mm = 110)

# ==============================================================================
# 10.5 FIGURE: METRIC-ASSOCIATION HEATMAP, BY GENUS AND OVERALL
# ==============================================================================
# A genus-resolved companion to the dataset-wide Akaike-weight importance.
# For every genus, and for the pooled dataset ("Overall"), the cell shows the
# Spearman correlation between each forest configuration metric and log(AOH).
# Spearman is used because n per genus is small and the metric-AOH relationship
# need not be linear within a genus. The colour encodes strength and direction
# (blue = negative, orange to red = positive, near-white = no association); the
# text annotation gives the correlation value.
#
# Important caveats, mirrored in the figure caption:
#   - This is descriptive, not inferential. Per-genus correlations with few
#     observations are unstable; cells with n < MIN_N_GENUS are greyed out as
#     "not estimable", and no p-values are reported per genus.
#   - The dataset-wide importance (sum of Akaike weights) remains the formal
#     importance metric (Section 10, Figure_PIII_03). The heatmap adds texture:
#     it shows whether a metric that matters overall matters consistently across
#     genera or is driven by a subset of them.
# ==============================================================================
message("\n[Visualization] Metric-association heatmap (by genus and overall)...")

# Long table of raw metric values with the response, one row per observation.
assoc_long <- data_p3 %>%
  select(GENUS, log_AOH, all_of(predictors_part3)) %>%
  pivot_longer(cols = all_of(predictors_part3),
               names_to = "Metric", values_to = "MetricValue")

# Per-genus Spearman correlations (only where n >= MIN_N_GENUS).
assoc_by_genus <- assoc_long %>%
  group_by(GENUS, Metric) %>%
  summarise(
    n   = sum(!is.na(MetricValue) & !is.na(log_AOH)),
    rho = if (n >= MIN_N_GENUS &&
              sd(MetricValue, na.rm = TRUE) > 0) {
      suppressWarnings(cor(MetricValue, log_AOH,
                           method = "spearman", use = "complete.obs"))
    } else NA_real_,
    .groups = "drop"
  ) %>%
  mutate(Group = as.character(GENUS)) %>%
  select(Group, Metric, n, rho)

# Overall (pooled) Spearman correlations.
assoc_overall <- assoc_long %>%
  group_by(Metric) %>%
  summarise(
    n   = sum(!is.na(MetricValue) & !is.na(log_AOH)),
    rho = suppressWarnings(cor(MetricValue, log_AOH,
                               method = "spearman", use = "complete.obs")),
    .groups = "drop"
  ) %>%
  mutate(Group = "Overall") %>%
  select(Group, Metric, n, rho)

assoc_df <- bind_rows(assoc_by_genus, assoc_overall)

# Export the underlying matrix for the Supplementary Material.
write.csv(assoc_df,
          file.path(output_path, "Table_S_PartIII_MetricAssociation_ByGenus.csv"),
          row.names = FALSE)

# Sample size per group, appended to each y-axis label for consistency with the
# effect heatmap (Section 10.6).
assoc_group_n <- assoc_df %>%
  group_by(Group) %>%
  summarise(n = max(n), .groups = "drop")
n_lookup_assoc  <- setNames(assoc_group_n$n, assoc_group_n$Group)
y_labeller_assoc <- function(x) paste0(x, " (n=", n_lookup_assoc[x], ")")

# Readable, single-line metric labels (the multi-line labels used elsewhere
# are too tall for a heatmap axis).
metric_labels_flat <- c(
  "PD_Forest"         = "Patch density (PD)",
  "Mean_FRAC_Forest"  = "Shape complexity (FRAC)",
  "Mean_ENN_Forest_m" = "Isolation (ENN)",
  "PROX_MN_Forest"    = "Proximity (PROX)",
  "PLAND_Forest"      = "Forest cover (PLAND)"
)

# Every predictor carried into the heatmaps must have a label. Stopping here is
# preferable to printing an axis that reads NA.
stopifnot(all(predictors_part3 %in% names(metric_labels_flat)))

assoc_plot_df <- assoc_df %>%
  mutate(
    Metric_lab = factor(metric_labels_flat[Metric],
                        levels = metric_labels_flat[predictors_part3]),
    # Keep "Overall" at the top; genera alphabetical below it.
    Group = factor(Group,
                   levels = c(sort(setdiff(unique(Group), "Overall")),
                              "Overall")),
    cell_label = ifelse(is.na(rho), "n.e.", sprintf("%.2f", rho)),
    # Text colour chosen per cell for legibility: dark on light/empty cells,
    # white on saturated (strongly coloured) cells.
    txt_col    = ifelse(is.na(rho) | abs(rho) < 0.55, "grey15", "white")
  )

plot_assoc <- ggplot(assoc_plot_df,
                     aes(x = Metric_lab, y = Group, fill = rho)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = cell_label, colour = txt_col),
            size = 3.1, fontface = "bold") +
  scale_colour_identity() +
  scale_fill_gradient2(low = heat_low, mid = heat_mid, high = heat_high,
                       midpoint = 0, limits = c(-1, 1),
                       na.value = "grey88",
                       name = expression("Spearman " * rho),
                       breaks = seq(-1, 1, 0.5)) +
  scale_y_discrete(labels = y_labeller_assoc) +
  labs(x = NULL, y = NULL) +
  coord_fixed() +
  theme_minimal(base_size = 11, base_family = "sans") +
  theme(plot.title          = element_text(face = "bold", size = 11, hjust = 0),
        plot.title.position = "plot",
        axis.text.x   = element_text(angle = 25, hjust = 1, size = 9),
        axis.text.y   = element_text(size = 9,
                                     face = ifelse(
                                       levels(assoc_plot_df$Group) == "Overall",
                                       "bold", "italic")),
        panel.grid    = element_blank(),
        legend.position = "right")

save_publication_plot(plot_assoc,
                      "Figure_PIII_07_Metric_Association_Heatmap.tiff",
                      width_mm = 168, height_mm = 150)

# ==============================================================================
# 10.6 FIGURE: Standardized EFFECT OF EACH PREDICTOR ON log(AOH),
#      BY GENUS AND FOR ALL SPECIES POOLED
# ==============================================================================
# This figure answers a direct, applied question: what is the effect of each
# forest configuration metric on the amount of habitat (log(AOH)), for all
# species together and for each genus separately?
#
# The estimate is deliberately the simplest one that is robust at this sample
# size: the Standardized slope (beta) of a univariate OLS, log(AOH) ~ metric_z,
# fitted within each genus and for all species pooled. The metric is on its
# global z-score, so beta is the change in log(AOH) per one standard deviation
# of the metric, on the same scale as the multivariable coefficients of
# Figure_PIII_02. The sign gives the direction, the magnitude gives the value of
# the variable, and the color encodes both (blue = negative, otange to red =
# positive).
#
# Why univariate, and not a model per genus. With 67 observations spread over
# the genera, each genus has only a handful of records. A four-predictor model,
# or an AICc model selection, fitted inside a single genus would be
# overparameterized and its coefficients unstable. A univariate standardized
# slope is the robust, honest per-genus estimate. It does not force a result:
# where a genus has too few observations (n < MIN_N_GENUS) or no variance in the
# metric, the cell is greyed out as "n.e." (not estimable).
#
# How to read it alongside the formal results.
#   - The dataset-wide, multivariable inference stays in Figure_PIII_02
#     (standardized coefficients of the final model) and Figure_PIII_03 (Akaike-
#     weight importance). Those use the full sample and remain the formal result.
#   - This heatmap is a descriptive companion that shows whether the pooled
#     effect is shared across genera or driven by a subset of them. Because the
#     per-genus slopes are univariate (one metric at a time) and the pooled
#     coefficient of Figure_PIII_02 is multivariable (all metrics jointly), a
#     univariate genus slope can differ in sign or magnitude from the
#     multivariable pooled effect, especially for PD_Forest and Mean_FRAC_Forest,
#     which are correlated. That difference is expected and informative.
# ==============================================================================
message("\n[Visualization] Standardised-effect heatmap (by genus and all species)...")

# Significance marker from a p-value (descriptive at the per-genus level).
sig_star <- function(p) {
  if (is.na(p))  return("")
  if (p < 0.001) return("***")
  if (p < 0.01)  return("**")
  if (p < 0.05)  return("*")
  ""
}

# Univariate standardised slope of log(AOH) on one metric, within a subset.
# metric_z is the GLOBAL z-score, so every beta is on a common "per global SD"
# scale and the genus rows and the pooled row are directly comparable.
univariate_slope <- function(df, metric) {
  z_var <- paste0(metric, "_z")
  x  <- df[[z_var]]
  y  <- df$log_AOH
  ok <- is.finite(x) & is.finite(y)
  x  <- x[ok]; y <- y[ok]
  n  <- length(y)
  if (n < MIN_N_GENUS || sd(x) == 0) {
    return(data.frame(n = n, beta = NA_real_, se = NA_real_, p = NA_real_))
  }
  fit <- lm(y ~ x)
  cf  <- summary(fit)$coefficients
  if (!"x" %in% rownames(cf)) {
    return(data.frame(n = n, beta = NA_real_, se = NA_real_, p = NA_real_))
  }
  data.frame(n    = n,
             beta = unname(cf["x", "Estimate"]),
             se   = unname(cf["x", "Std. Error"]),
             p    = unname(cf["x", "Pr(>|t|)"]))
}

genera_vec <- sort(as.character(unique(data_p3$GENUS)))
groups_vec <- c(genera_vec, "All species")

effect_rows <- list()
for (g in groups_vec) {
  sub <- if (g == "All species") data_p3 else dplyr::filter(data_p3, as.character(GENUS) == g)
  for (m in predictors_part3) {
    s <- univariate_slope(sub, m)
    effect_rows[[length(effect_rows) + 1L]] <- data.frame(
      Group  = g,
      Metric = m,
      n      = s$n,
      beta   = round(s$beta, 4),
      se     = round(s$se, 4),
      p      = round(s$p, 4),
      stringsAsFactors = FALSE
    )
  }
}
effect_df <- dplyr::bind_rows(effect_rows)

# Export the matrix that underlies the figure.
write.csv(effect_df,
          file.path(output_path, "Table_S_PartIII_EffectByGenus.csv"),
          row.names = FALSE)

# --- 10.6b Genus-specific slopes with the genus as a moderator ----------------
# The block above fits one model inside each genus, so a genus with four units
# estimates its own residual variance from four units. This block asks the same
# question in the opposite way: one model per metric over all 67 units, with the
# genus entered as a moderator, so every genus slope is read off a single fitted
# surface and shares one residual variance. Genera with few units therefore
# borrow precision from the rest, and the two tables bracket the per-genus
# result rather than duplicating it. Where they disagree, the disagreement is
# about how much a small genus is allowed to speak for itself, and both are
# reported for that reason.
#
# The response is kept on the observed scale and the family is Gamma with a
# logarithmic link, which is the family used elsewhere in this script, so the
# coefficients are on the same multiplicative scale as the reported model. Each
# genus slope is the linear contrast that adds the reference slope to the
# interaction term of that genus. The standard error comes from the full
# covariance matrix of the fit, not from the diagonal alone, because the two
# terms are correlated. Intervals use the t quantile on the residual degrees of
# freedom, which is the convention this script applies to every Gamma fit.
#
# This block replaces the marginal-effects script that was removed. It uses only
# base R, so no additional package is attached to obtain it.

slope_by_genus_moderator <- function(dat, metric) {
  x <- as.numeric(scale(dat[[metric]]))
  if (all(is.na(x)) || stats::sd(x, na.rm = TRUE) == 0) return(NULL)
  d_fit <- data.frame(y = dat[["AOH_ha"]], x = x,
                      GENUS = factor(as.character(dat$GENUS)))
  d_fit <- d_fit[stats::complete.cases(d_fit), ]
  fit <- stats::glm(y ~ x * GENUS, family = stats::Gamma(link = "log"),
                    data = d_fit)
  bt <- stats::coef(fit)
  V  <- stats::vcov(fit)
  tq <- stats::qt(0.975, stats::df.residual(fit))
  levs <- levels(d_fit$GENUS)
  out <- lapply(levs, function(g) {
    cv <- rep(0, length(bt)); names(cv) <- names(bt)
    cv["x"] <- 1
    inter <- paste0("x:GENUS", g)
    if (inter %in% names(bt)) cv[inter] <- 1
    b  <- sum(cv * bt)
    se <- sqrt(as.numeric(t(cv) %*% V %*% cv))
    data.frame(
      Metric   = metric,
      Genus    = g,
      N_units  = sum(d_fit$GENUS == g),
      Slope    = round(b, 4),
      CI_lower = round(b - tq * se, 4),
      CI_upper = round(b + tq * se, 4),
      stringsAsFactors = FALSE
    )
  })
  dplyr::bind_rows(out)
}

moderator_df <- dplyr::bind_rows(
  lapply(predictors_part3, function(m) slope_by_genus_moderator(data_p3, m))
)
moderator_df$Significant <- ifelse(
  moderator_df$CI_lower > 0 | moderator_df$CI_upper < 0, "yes", "no")

# The count that the supplementary caption states, computed rather than typed.
n_excluding_zero <- sum(moderator_df$Significant == "yes")
message(sprintf(
  "Moderator slopes: %d of %d intervals exclude zero.",
  n_excluding_zero, nrow(moderator_df)))

write.csv(moderator_df,
          file.path(output_path,
                    "Table_S_PartIII_SlopeByGenus_Moderator.csv"),
          row.names = FALSE)

# --- Assemble the single-panel heatmap ----------------------------------------
# "All species" is the last factor level, so it sits at the top of the y axis
# and is set in bold; genera are alphabetical below it.
group_levels <- c(genera_vec, "All species")

# Sample size per group, appended to each y-axis label so the reader can weigh
# each row by its support (small genera carry less weight).
group_n   <- effect_df %>%
  dplyr::group_by(Group) %>%
  dplyr::summarise(n = max(n), .groups = "drop")
n_lookup  <- setNames(group_n$n, group_n$Group)
y_labeller <- function(x) paste0(x, " (n=", n_lookup[x], ")")

# Symmetric colour limit for the diverging effect scale.
beta_lim <- max(abs(effect_df$beta), na.rm = TRUE)
beta_lim <- ifelse(is.finite(beta_lim), ceiling(beta_lim * 10) / 10, 1)

effect_plot_df <- effect_df %>%
  dplyr::mutate(
    Metric_lab = factor(metric_labels_flat[Metric],
                        levels = metric_labels_flat[predictors_part3]),
    Group      = factor(Group, levels = group_levels),
    cell_label = ifelse(is.na(beta), "n.e.",
                        paste0(sprintf("%.2f", beta),
                               vapply(p, sig_star, character(1)))),
    # Dark text on light cells, white text on strongly coloured cells.
    txt_col    = ifelse(!is.na(beta) & abs(beta) > 0.6 * beta_lim,
                        "white", "grey15")
  )

plot_effect_heatmap <- ggplot(effect_plot_df,
                              aes(x = Metric_lab, y = Group, fill = beta)) +
  geom_tile(color = "white", linewidth = 0.8) +
  geom_text(aes(label = cell_label, colour = txt_col),
            size = 3.1, fontface = "bold") +
  scale_colour_identity() +
  scale_fill_gradient2(low = heat_low, mid = heat_mid, high = heat_high,
                       midpoint = 0, limits = c(-beta_lim, beta_lim),
                       na.value = "grey88",
                       name = expression("Standardised " * beta),
                       breaks = scales::pretty_breaks(n = 5)) +
  scale_y_discrete(labels = y_labeller) +
  labs(x = NULL, y = NULL) +
  coord_fixed() +
  theme_minimal(base_size = 11, base_family = "sans") +
  theme(plot.title          = element_text(face = "bold", size = 11, hjust = 0),
        plot.title.position = "plot",
        axis.text.x   = element_text(angle = 25, hjust = 1, size = 9),
        axis.text.y   = element_text(size = 9,
                                     face = ifelse(group_levels == "All species",
                                                   "bold", "italic")),
        panel.grid    = element_blank(),
        legend.position = "right")

save_publication_plot(plot_effect_heatmap,
                      "Figure_PIII_08_Effect_Heatmap.tiff",
                      width_mm = 168, height_mm = 150)

# ==============================================================================
# 11. FIGURE: PARTIAL EFFECTS (predicted response curves)
# ==============================================================================
# Predicted response for each predictor while the others are held at their
# mean value (0 on the z-score scale). Predictions are back-transformed to the
# original units of AOH (ha) for biological interpretation.
message("\n[Visualization] Partial effect curves...")

# Choose the model used for prediction.
# When the competitive set contains a single model, that model is used directly
# for both inference and prediction. When model averaging is invoked, the
# coefficient forest plot already reports the averaged coefficients (Section
# 9), and predictions for the partial-effect curves are generated from the
# highest-weight competitive model. This is a deliberate, conservative choice:
# MuMIn::predict.averaging is sensitive to mismatches between dredge term
# names and the new data frame, and using the best-supported single model for
# prediction sidesteps that issue while preserving the multimodel character
# of the inference itself (Burnham & Anderson 2002).
predict_model <- if (averaged) competitive[[1]] else final_model
is_glm_pred   <- inherits(predict_model, "glm") &&
  family(predict_model)$family == "Gamma"

# Build a long data frame of predictions on a grid for each predictor.
# Both branches predict on the link scale (log), build the 95% confidence
# interval there, and only then back-transform via exp(). This is the correct
# protocol for log-scale models and avoids the asymmetry pitfalls of
# back-transforming on the response scale.
make_partial_df <- function(varname, model, data, raw_data) {
  z_var   <- paste0(varname, "_z")
  grid_z  <- seq(min(data[[z_var]], na.rm = TRUE),
                 max(data[[z_var]], na.rm = TRUE),
                 length.out = 100)
  newdata <- data %>%
    summarise(across(ends_with("_z"), ~0)) %>%
    slice(rep(1, length(grid_z)))
  newdata[[z_var]] <- grid_z
  
  # For OLS on log(AOH): predict on response scale returns log-AOH predictions.
  # For Gamma GLM with log link: predict on link scale returns log-mu directly.
  # ADDED 30 Jul 2026. predict.merMod has no se.fit argument, and its model
  # frame requires the grouping variables, which the new data built above does
  # not carry because it keeps only the "_z" columns. That is the source of the
  # error "object 'GENUS' not found". For a mixed model the population-level
  # (marginal) curve is therefore built directly from the fixed-effect design
  # matrix and its covariance, which is what a partial-effect plot should show:
  # the effect of the predictor with the random intercept set to its mean of
  # zero, not the fit conditional on any one genus.
  if (inherits(model, "merMod")) {
    for (g in names(lme4::getME(model, "flist"))) {
      newdata[[g]] <- data[[g]][1]          # placeholder level, unused below
    }
    fe_terms <- stats::delete.response(
      stats::terms(lme4::nobars(stats::formula(model))))
    X        <- stats::model.matrix(fe_terms, data = newdata)
    beta     <- lme4::fixef(model)
    X        <- X[, names(beta), drop = FALSE]
    V        <- as.matrix(stats::vcov(model))
    link_fit <- as.vector(X %*% beta)
    se_fit   <- sqrt(pmax(rowSums((X %*% V) * X), 0))
  } else if (is_glm_pred) {
    pred     <- predict(model, newdata = newdata, type = "link", se.fit = TRUE)
    link_fit <- pred$fit
    se_fit   <- pred$se.fit
  } else {
    pred     <- predict(model, newdata = newdata, se.fit = TRUE)
    link_fit <- pred$fit
    se_fit   <- pred$se.fit
  }
  link_lo   <- link_fit - 1.96 * se_fit
  link_up   <- link_fit + 1.96 * se_fit
  
  data.frame(
    Predictor = varname,
    x_real    = grid_z * sd(raw_data[[varname]], na.rm = TRUE) +
      mean(raw_data[[varname]], na.rm = TRUE),
    fit       = exp(link_fit),
    ci_low    = exp(link_lo),
    ci_up     = exp(link_up)
  )
}

partial_df <- bind_rows(lapply(predictors_part3, make_partial_df,
                               model    = predict_model,
                               data     = data_p3,
                               raw_data = data_p3)) %>%
  mutate(Predictor_label = label_predictors(Predictor))

raw_points_df <- data_p3 %>%
  select(AOH_ha, all_of(predictors_part3)) %>%
  pivot_longer(cols = all_of(predictors_part3),
               names_to = "Predictor", values_to = "x_real") %>%
  mutate(Predictor_label = label_predictors(Predictor))

plot_partial <- ggplot() +
  geom_point(data = raw_points_df,
             aes(x = x_real, y = AOH_ha),
             alpha = 0.35, size = 1.6, color = "grey45") +
  geom_ribbon(data = partial_df,
              aes(x = x_real, ymin = ci_low, ymax = ci_up,
                  fill = Predictor), alpha = 0.25, show.legend = FALSE) +
  geom_line(data = partial_df,
            aes(x = x_real, y = fit, color = Predictor),
            linewidth = 1.1, show.legend = FALSE) +
  facet_wrap(~ Predictor_label, scales = "free_x", ncol = 3) +
  scale_color_manual(values = colors_predictors) +
  scale_fill_manual(values  = colors_predictors) +
  scale_y_log10(labels = label_comma()) +
  labs(x = NULL,
       y = "Predicted OHE (ha, log scale)") +
  theme_classic(base_size = 10, base_family = "sans") +
  theme(plot.title = element_text(face = "bold", size = 11, hjust = 0.5),
        strip.background = element_blank(),
        strip.text = element_text(face = "bold", size = 9))

save_publication_plot(plot_partial,
                      "Figure_PIII_04_Partial_Effects.tiff",
                      width_mm = 174, height_mm = 130)

# ==============================================================================
# 12. FIGURE: OBSERVED vs PREDICTED
# ==============================================================================
# Diagnostic plot used to verify that the final model captures the gradient of
# AOH, not only the central tendency.
message("\n[Visualization] Observed vs predicted...")

obs_pred_df <- data_p3 %>%
  mutate(predicted = if (is_glm_pred) {
    predict(predict_model, type = "response")
  } else if (inherits(predict_model, "merMod")) {
    # re.form = NA gives the population-level prediction, consistent with the
    # marginal partial-effect curves plotted above.
    exp(predict(predict_model, re.form = NA))
  } else {
    exp(predict(predict_model, type = "response"))
  })

r_squared_obs_pred <- cor(obs_pred_df$AOH_ha, obs_pred_df$predicted)^2

# ADDED 23 Aug 2026. The caption of Figure S5 quotes the value annotated in the
# panel and its counterpart on the logarithmic scale. Both are exported here so
# that the caption is traceable to an output file rather than read off the image.
r_squared_obs_pred_log <- cor(log(obs_pred_df$AOH_ha),
                              log(obs_pred_df$predicted))^2
write.csv(
  data.frame(
    Quantity = c("Squared correlation, observed against predicted, hectare scale",
                 "Squared correlation, observed against predicted, logarithmic scale",
                 "Model used for the figure",
                 "Largest sampling unit (ha)"),
    Value    = c(round(r_squared_obs_pred, 4),
                 round(r_squared_obs_pred_log, 4),
                 NA, round(max(obs_pred_df$AOH_ha), 2)),
    Note     = c("annotated in the panel", "",
                 model_family_label, ""),
    stringsAsFactors = FALSE
  ),
  file.path(output_path, "Table_S_PartIII_ObsPred_R2.csv"),
  row.names = FALSE, na = "")
message(sprintf("   R2 obs-pred: %.4f on the hectare scale, %.4f on the log scale",
                r_squared_obs_pred, r_squared_obs_pred_log))

plot_obs_pred <- ggplot(obs_pred_df, aes(x = predicted, y = AOH_ha)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              color = "grey40", linewidth = 0.6) +
  geom_point(alpha = 0.55, size = 2, color = "#0072B2") +
  geom_smooth(method = "lm", se = TRUE,
              color = "#D55E00", fill = "#FED9B7", linewidth = 0.7) +
  scale_x_log10(labels = label_comma()) +
  scale_y_log10(labels = label_comma()) +
  annotate("text", x = min(obs_pred_df$predicted), y = max(obs_pred_df$AOH_ha),
           label = sprintf("R^2 == %.3f", r_squared_obs_pred),
           parse = TRUE, hjust = 0, vjust = 1, size = 4) +
  labs(x = "Predicted OHE (ha)", y = "Observed OHE (ha)") +
  theme_classic(base_size = 11, base_family = "sans") +
  theme(plot.title = element_text(face = "bold", size = 11, hjust = 0.5))

save_publication_plot(plot_obs_pred,
                      "Figure_PIII_05_Observed_vs_Predicted.tiff",
                      width_mm = 130, height_mm = 130)

# ==============================================================================
# 13. FIGURE: RESIDUAL DIAGNOSTIC PANEL (visual)
# ==============================================================================
# Four-panel diagnostic figure for the OLS or GLM. Complements the formal
# Shapiro / Breusch-Pagan / DHARMa tests reported in Section 6.
message("\n[Visualization] Residual diagnostic panel...")

tiff(file.path(output_path, "Figure_PIII_06_Residual_Diagnostics.tiff"),
     width = 174 * 600 / 25.4, height = 174 * 600 / 25.4, units = "px",
     res = 600, compression = "lzw")
par(mfrow = c(2, 2), mar = c(4, 4, 2.5, 1), oma = c(0, 0, 1.5, 0),
    family = "sans", cex.lab = 0.95, cex.main = 1.0)
# ADDED 30 Jul 2026. plot.merMod is lattice-based and does not accept the
# which= argument of plot.lm, so the four panels are drawn by hand when the
# retained model is a mixed one. Panels 1 to 3 mirror plot.lm; the fourth
# replaces residuals-versus-leverage, which is not defined in the same way for
# a mixed model, by the distribution of residuals across the grouping factor,
# which is the more informative check at this level.
if (inherits(global_model, "merMod")) {
  res_g   <- stats::residuals(global_model)
  fit_g   <- stats::fitted(global_model)
  sres_g  <- res_g / stats::sd(res_g)
  grp_g   <- lme4::getME(global_model, "flist")[[1]]

  plot(fit_g, res_g, xlab = "Fitted values", ylab = "Residuals",
       main = "Residuals vs fitted", pch = 19, col = "#00000088")
  abline(h = 0, lty = 2, col = "grey40")
  lines(stats::lowess(fit_g, res_g), col = "#D55E00", lwd = 2)

  stats::qqnorm(sres_g, main = "Normal Q-Q", pch = 19, col = "#00000088")
  stats::qqline(sres_g, lty = 2, col = "grey40")

  plot(fit_g, sqrt(abs(sres_g)), xlab = "Fitted values",
       ylab = expression(sqrt(abs("standardised residuals"))),
       main = "Scale-location", pch = 19, col = "#00000088")
  lines(stats::lowess(fit_g, sqrt(abs(sres_g))), col = "#D55E00", lwd = 2)

  boxplot(res_g ~ grp_g, xlab = "", ylab = "Residuals",
          main = "Residuals by grouping level", las = 2,
          cex.axis = 0.7, col = "grey90")
  abline(h = 0, lty = 2, col = "grey40")
} else {
  # sub.caption = "" suppresses the deparsed call that plot.lm writes into the
  # outer margin, which overprinted the title added by mtext() below.
  plot(global_model, which = c(1, 2, 3, 5), sub.caption = "")
}
mtext(sprintf("Residual diagnostics: %s", model_family_label),
      outer = TRUE, cex = 1.0, font = 2)
invisible(dev.off())
message("   Residual diagnostic panel exported.")

# ==============================================================================
# 14. CONSOLE SUMMARY
# ==============================================================================
message("\n============================================================")
message("  PART III — Summary of results")
message("============================================================")
message(sprintf("  Observations / SUs            : %d / %d", n_obs, n_su))
message(sprintf("  Inferential framework         : %s", model_family_label))
message(sprintf("  Global-model R^2 / DevExpl    : %.3f",
                glob_summary$R2_or_DevExpl))
if (exists("icc_su")) {
  message(sprintf("  ICC of the genus              : %.3f", icc_su))
}
message(sprintf("  Competitive set (Delta < 2)   : %d model(s)", n_competitive))
message(sprintf("  Inference strategy            : %s", inf_label))
message("\n  Predictor importance (Akaike weights):")
for (i in seq_len(nrow(df_importance))) {
  message(sprintf("    - %-22s : %.3f",
                  df_importance$Predictor[i], df_importance$Importance[i]))
}
message("============================================================\n")

# ==============================================================================
# 15. EXPORT MANIFEST
# ==============================================================================
message("[Export] Files written by Part III:")
message("  Tables:")
message("    - Table_S_PartIII_Averaging_FullVsConditional.csv")
message("    - Table_S_PartIII_CompetitiveModels.csv")
message("    - Table_S_PartIII_ConfirmatoryModel.csv")
message("    - Table_S_PartIII_DependenceDiagnostic.csv")
message("    - Table_S_PartIII_DescriptiveStatistics.csv")
message("    - Table_S_PartIII_Diagnostics.csv")
message("    - Table_S_PartIII_Dredge_FullTable.csv")
message("    - Table_S_PartIII_EffectByGenus.csv")
message("    - Table_S_PartIII_FRAC_Area_Correlation.csv")
message("    - Table_S_PartIII_FinalCoefficients.csv")
message("    - Table_S_PartIII_GlobalModelSummary.csv")
message("    - Table_S_PartIII_MetricAssociation_ByGenus.csv")
message("    - Table_S_PartIII_ModelComparison.csv")
message("    - Table_S_PartIII_ObsPred_R2.csv")
message("    - Table_S_PartIII_Sensitivity_MatrixModel.csv")
message("    - Table_S_PartIII_SlopeByGenus_Moderator.csv")
message("    - Table_S_PartIII_SpatialSensitivity.csv")
message("    - Table_S_PartIII_StabilityCheck.csv")
message("    - Table_S_PartIII_VIF.csv")
message("    - Table_S_PartIII_VariableImportance.csv")
message("  Figures (TIFF, 600 dpi):")
message("    - Figure_PIII_01_Response_Transformation.tiff")
message("    - Figure_PIII_02_Coefficient_ForestPlot.tiff")
message("    - Figure_PIII_03_Variable_Importance.tiff")
message("    - Figure_PIII_04_Partial_Effects.tiff")
message("    - Figure_PIII_05_Observed_vs_Predicted.tiff")
message("    - Figure_PIII_06_Residual_Diagnostics.tiff")
message("    - Figure_PIII_07_Metric_Association_Heatmap.tiff")
message("    - Figure_PIII_08_Effect_Heatmap.tiff")
message("    - session_info_part03.txt")

# ==============================================================================
# 16. WHERE THE METHODOLOGICAL REASONING IS RECORDED
# ==============================================================================
# This file carried, until 24 August 2026, a block of interpretive notes written
# for an earlier version of the study. It was removed because it had ceased to
# describe this script: it called the response AOH, which the manuscript now
# calls the occupied habitat envelope; it stated that the 67 observations fall
# across 12 sampling units and that the mixed specification uses a random
# intercept for the sampling unit, whereas there are 67 units and the random
# intercept is for the genus; and it quoted variance fractions, correlations and
# a mapping accuracy that the chain has since recomputed. Notes that contradict
# the code they sit in are worse than no notes.
#
# The reasoning now lives where it can be checked against a number:
#
#   Manuscript, Section 2.5.1   why three metrics were excluded, and why a
#                               variance inflation factor does not detect the
#                               dependence that excluded them
#   Manuscript, Section 2.5.5   the six candidate specifications and the order
#                               in which the three criteria are applied
#   Manuscript, Section 4.1     what the sampling window does to a normalized
#                               metric, and which coefficients survive it
#   docs/TABLE_MAP.md           which script writes each table and figure of the
#                               manuscript and of the supplement
#   docs/SCRIPTS.md             what each script does, in order
#   docs/Registro_de_Correcoes_22.08.md and _23.08.md
#                               every defect found and corrected, with the value
#                               before and after
#
# The archived text of the old block is kept outside the repository, with the
# revision record of 24 August 2026.
# ==============================================================================

# ==============================================================================
# 17. SESSION INFO
# ==============================================================================
message("\n[Session] Package versions and system information:")
print(sessionInfo())
writeLines(capture.output(sessionInfo()),
           file.path(output_path, "session_info_part03.txt"))

# ==============================================================================
# REFERENCES (key methodological literature)
#
# Bolker, B.M., Brooks, M.E., Clark, C.J., Geange, S.W., Poulsen, J.R., Stevens,
#   M.H.H., White, J.S.S. (2009). Generalized linear mixed models: a practical
#   guide for ecology and evolution. Trends in Ecology & Evolution, 24, 127-135.
# Borcard, D., Gillet, F., Legendre, P. (2018). Numerical Ecology with R, 2nd
#   ed. Springer.
# Brooks, T.M., Pimm, S.L., Akcakaya, H.R., Buchanan, G.M., Butchart, S.H.M.,
#   Foden, W., Hilton-Taylor, C., Hoffmann, M., Jenkins, C.N., Joppa, L.,
#   Li, B.V., Menon, V., Ocampo-Penuela, N., Rondinini, C. (2019). Measuring
#   terrestrial Area of Habitat (AOH) and its utility for the IUCN Red List.
#   Trends in Ecology & Evolution, 34, 977-986.
# Burnham, K.P., Anderson, D.R. (2002). Model Selection and Multimodel
#   Inference: A Practical Information-Theoretic Approach, 2nd ed. Springer.
# Dormann, C.F. et al. (2013). Collinearity: a review of methods to deal
#   with it and a simulation study evaluating their performance. Ecography,
#   36, 27-46.
# Grueber, C.E., Nakagawa, S., Laws, R.J., Jamieson, I.G. (2011). Multimodel
#   inference in ecology and evolution: challenges and solutions. Journal of
#   Evolutionary Biology, 24, 699-711.
# Nakagawa, S., Schielzeth, H. (2013). A general and simple method for
#   obtaining R2 from generalized linear mixed-effects models. Methods in
#   Ecology and Evolution, 4, 133-142.
# Zuur, A.F., Ieno, E.N., Walker, N., Saveliev, A.A., Smith, G.M. (2009). Mixed
#   Effects Models and Extensions in Ecology with R. Springer.
# Zuur, A.F., Ieno, E.N., Elphick, C.S. (2010). A protocol for data exploration
#   to avoid common statistical problems. Methods in Ecology and Evolution,
#   1, 3-14.
# ==============================================================================