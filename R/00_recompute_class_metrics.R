# ==============================================================================
# 00  Recomputation of the class-level landscape metrics
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
#   Recompute every class-level metric for all 67 sampling units directly from
#   the per-unit rasters, and compare the result with the values already stored.
#   Two units reached the consolidated dataset with empty cells; without this
#   step complete.cases() drops them silently further downstream. Nothing is
#   imputed: every filled value is measured again from the same raster.
#
# INPUT
#   67 per-unit rasters; FRAGSTATS patch tables for the proximity index;
#   Dados/Processados/Data_Raw_FINAL.csv for the agreement check
#
# OUTPUT
#   Dados/Processados/Tabela_Metricas_Recomputada.csv
#   Outputs/Manuscrito/Part0/TableS_Recompute_Agreement.csv
#   Outputs/Manuscrito/Part0/TableS_Repaired_Cells.csv
#
# RUNTIME
#   Minutes. Each raster is cropped with terra::trim() before any metric is
#   computed, because more than 99 per cent of each file is NoData padding, and
#   a one-cell margin is added back with terra::extend() so that edge density
#   cannot react to the crop. Per-unit results are cached, so an interrupted run
#   resumes.
# ==============================================================================


suppressPackageStartupMessages({
  library(here)
  library(terra)
  library(landscapemetrics)
  library(dplyr)
  library(stringr)
  library(purrr)
  library(readr)
  library(tidyr)
})

GLOBAL_SEED <- 123
set.seed(GLOBAL_SEED)

# ------------------------------------------------------------------------------
# RUN CONTROL
# ------------------------------------------------------------------------------
# RECOMPUTE_PATCH_MEANS. FALSE computes only the metrics needed to repair the
# dataset and to check agreement on the repaired columns: CA, PLAND, NP, PD, ED
# and LPI. TRUE adds AREA_MN, FRAC_MN and ENN_MN, which are not repaired because
# they were never empty, and which cost most of the run time. ENN_MN alone is
# quadratic in the number of patches.
RECOMPUTE_PATCH_MEANS <- FALSE

# USE_CACHE. TRUE stores each unit's result under cache_path and skips units
# already computed. Delete the directory to force a clean recomputation.
USE_CACHE <- TRUE

raster_path    <- here("Dados", "Rasters", "UA_RASTER")
fragstats_path <- here("Dados", "FRAGSTATS_RESULT")
out_path       <- here("Dados", "Processados")
report_path    <- here("Outputs", "Manuscrito", "Part0")
cache_path     <- here("Dados", "Processados", "_part0_cache")

if (!dir.exists(report_path)) dir.create(report_path, recursive = TRUE)
if (USE_CACHE && !dir.exists(cache_path)) dir.create(cache_path, recursive = TRUE)

N_UNITS_EXPECTED <- 67

# Class codes of the MapBiomas aggregation used throughout the study.
# Farming is the union of cls_3 and cls_4, so cls_4 is reclassified into cls_3
# BEFORE the metrics are computed. This is what makes NP, PD and LPI correct for
# that class: adjacent pasture and agriculture pixels form a single patch.
CLASS_CODES    <- c(Forest = 1, Herbaceous = 2, Agropecuaria = 3,
                    Water = 5, Non_Vegetated = 6)
FARMING_FROM   <- 4
FARMING_TO     <- 3


# ==============================================================================
# 2. TAXON NAME Normalization
# ==============================================================================
# Every orthographic deviation found in the project is listed here, so that the
# correction is explicit and auditable rather than silent. Three of them are
# genuine typographical errors and one is a nomenclatural correction.
#
#   aguriba        raster ud6_aguriba1.tif, missing the second 'a' of guariba.
#                  This is the file that dropped A_guariba UA 6.
#   baracnoides    missing the 'h' of arachnoides, in the intermediate table.
#                  This is the label that dropped B_arachnoides UA 1.
#   concolor       Puma rasters ud10 to ud12 use 'concolor' while ud1 to ud9 use
#                  'pconcolor'. Both denote Puma concolor.
#   B_hypoxantus   The valid spelling of the taxon is Brachyteles hypoxanthus
#                  (Kuhl, 1820). The tabular data carried B_hypoxantus, without
#                  the 'h', while the shapefile carried the correct spelling in
#                  SP_ID. The compound key therefore never matched and the two
#                  units of the species reached Part II without centroid
#                  coordinates. The valid spelling is adopted here and in Parts
#                  I, II and III.
#
# canonical_taxon() is reproduced in Parts I, II and III. If a name is added here
# it must be added there as well. Until 24 August 2026 this copy returned NA for
# an unrecognised name while the other three stopped, so a name absent from the
# list propagated silently into the join instead of halting the run.
# ==============================================================================
canonical_taxon <- function(x) {
  key <- x %>%
    as.character() %>%
    str_to_lower() %>%
    str_replace_all("[^a-z]", "")   # drops underscores, spaces, digits, dots

  dplyr::case_when(
    key %in% c("aguariba", "aguriba", "alouattaguariba")        ~ "A_guariba",
    key %in% c("barach", "baracnoides", "barachnoides",
               "brachytelesarachnoides")                        ~ "B_arachnoides",
    key %in% c("bcrinitus", "btorquatus", "bradypuscrinitus")    ~ "B_crinitus",
    key %in% c("bhypox", "bhypoxantus", "bhypoxanthus",
               "brachyteleshypoxanthus")                        ~ "B_hypoxanthus",
    key %in% c("lwiedii", "leoparduswiedii")                     ~ "L_wiedii",
    key %in% c("mazama", "mazamasp")                             ~ "Mazama",
    key %in% c("mtrida", "mtridactyla",
               "myrmecophagatridactyla")                         ~ "M_tridactyla",
    key %in% c("pconcolor", "concolor", "pumaconcolor")          ~ "P_concolor",
    key %in% c("tpecari", "tayassupecari")                       ~ "T_pecari",
    TRUE ~ NA_character_
  ) -> out

  if (any(is.na(out) & !is.na(x))) {
    stop("Unrecognised taxon name(s): ",
         paste(unique(x[is.na(out) & !is.na(x)]), collapse = ", "),
         ". Add them to canonical_taxon() in Parts 0, I, II and III before ",
         "rerunning.", call. = FALSE)
  }
  out
}


# ==============================================================================
# 3. READING THE RASTERS                                               [1/5]
# ==============================================================================
message("\n[1/5] Indexing the per-unit rasters...")

files <- list.files(raster_path, pattern = "\\.tif$", full.names = TRUE)
message(sprintf("   %d rasters found in %s", length(files), raster_path))

if (length(files) != N_UNITS_EXPECTED) {
  stop(sprintf(paste0("Expected %d rasters, found %d. The recomputation must ",
                      "cover every sampling unit before it is used to repair ",
                      "the dataset."),
               N_UNITS_EXPECTED, length(files)))
}

keys <- tibble(
  path    = files,
  file    = basename(files),
  # File pattern: ud<UA>_<taxon><n>.tif. The unit is the run of digits after 'ud'.
  UA_ID   = as.integer(str_match(basename(files), "^ud(\\d+)_")[, 2]),
  raw_tax = str_match(basename(files), "^ud\\d+_([a-z]+)\\d*\\.tif$")[, 2]
) %>%
  mutate(SPECIES = canonical_taxon(raw_tax))

if (any(is.na(keys$SPECIES)) || any(is.na(keys$UA_ID))) {
  print(keys %>% filter(is.na(SPECIES) | is.na(UA_ID)), n = Inf)
  stop("Unparsed file names. Add the pattern to canonical_taxon() before rerunning.")
}
if (anyDuplicated(keys[, c("SPECIES", "UA_ID")]) > 0) {
  print(keys[duplicated(keys[, c("SPECIES", "UA_ID")]), ], n = Inf)
  stop("Duplicated (SPECIES, UA_ID) after normalisation.")
}

message(sprintf("   %d unique (taxon, unit) pairs across %d taxa.",
                nrow(keys), n_distinct(keys$SPECIES)))
print(keys %>% count(SPECIES, name = "n_units"), n = Inf)


# ==============================================================================
# 4. CLASS-LEVEL METRICS FROM THE RASTERS                              [2/5]
# ==============================================================================
# directions = 8 and the default count_boundary = FALSE of lsm_c_ed reproduce
# the settings of Calculos_Landscape_Metrics.R. Do not change them here without
# recomputing the whole column for all 67 units.
# ==============================================================================
message("\n[2/5] Computing class-level metrics with landscapemetrics...")

# Metric set. The core is what the repair needs and what the agreement check must
# cover, because ED, LPI, NP and PD are precisely the columns that are empty in
# the two units. The three per-patch means are optional; see RECOMPUTE_PATCH_MEANS.
WHAT_CORE <- c("lsm_c_ca", "lsm_c_pland", "lsm_c_np", "lsm_c_pd",
               "lsm_c_ed", "lsm_c_lpi")
WHAT_SLOW <- c("lsm_c_area_mn", "lsm_c_frac_mn", "lsm_c_enn_mn")
WHAT_USE  <- if (RECOMPUTE_PATCH_MEANS) c(WHAT_CORE, WHAT_SLOW) else WHAT_CORE

message(sprintf("   Metrics requested: %s", paste(WHAT_USE, collapse = ", ")))
if (!RECOMPUTE_PATCH_MEANS) {
  message("   AREA_MN, FRAC_MN and ENN_MN skipped (RECOMPUTE_PATCH_MEANS = FALSE).")
  message("   They are not repaired by this script and ENN_MN dominates the run time.")
}

unit_metrics <- function(path, sp, ua) {
  r <- terra::rast(path)
  n_before <- terra::ncell(r)

  # Merge cls_4 into cls_3 so that the farming class is a single class and its
  # patch counts reflect the merged geometry.
  r <- terra::classify(r, cbind(FARMING_FROM, FARMING_TO))

  # Crop to the bounding box of the valid pixels. The per-unit rasters cover the
  # full extent of every unit of the taxon, so typically more than 99 per cent of
  # the file is NoData padding. See the header for the measured reduction.
  r <- terra::trim(r)

  # Put a one-cell NoData margin back, so that every valid pixel keeps exactly
  # the neighbourhood it had in the original file and Edge Density cannot react
  # to the crop. Section 6 verifies this against the stored values.
  r <- terra::extend(r, 1)

  n_after <- terra::ncell(r)

  # freq() counts valid cells without materialising the raster as an R vector.
  fr      <- terra::freq(r)
  area_ha <- sum(fr$count) * prod(terra::res(r)) / 10000

  res <- landscapemetrics::calculate_lsm(
    r,
    level      = "class",
    what       = WHAT_USE,
    directions = 8,
    verbose    = FALSE
  )

  attr(res, "n_before") <- n_before
  attr(res, "n_after")  <- n_after

  res %>%
    mutate(SPECIES = sp, UA_ID = ua, UA_Area_ha_recomp = area_ha) %>%
    select(SPECIES, UA_ID, UA_Area_ha_recomp, class, metric, value)
}

# Cache key. The flag is part of the name so that switching RECOMPUTE_PATCH_MEANS
# does not silently reuse a result computed with the other metric set.
cache_file <- function(sp, ua) {
  file.path(cache_path,
            sprintf("%s_UA%02d_%s.rds", sp, ua,
                    if (RECOMPUTE_PATCH_MEANS) "full" else "core"))
}

n_units  <- nrow(keys)
t_start  <- Sys.time()
out_list <- vector("list", n_units)
n_cached <- 0L

for (i in seq_len(n_units)) {
  sp <- keys$SPECIES[i]; ua <- keys$UA_ID[i]; pth <- keys$path[i]
  cf <- cache_file(sp, ua)

  if (USE_CACHE && file.exists(cf)) {
    out_list[[i]] <- readRDS(cf)
    n_cached <- n_cached + 1L
    message(sprintf("   [%2d/%2d] %-16s UA %-3d  cached", i, n_units, sp, ua))
    next
  }

  t0  <- Sys.time()
  res <- unit_metrics(pth, sp, ua)
  dt  <- as.numeric(difftime(Sys.time(), t0, units = "secs"))

  if (USE_CACHE) saveRDS(res, cf)
  out_list[[i]] <- res

  done      <- i - n_cached
  elapsed   <- as.numeric(difftime(Sys.time(), t_start, units = "secs"))
  remaining <- if (done > 0) (elapsed / done) * (n_units - i) else NA_real_

  message(sprintf(
    "   [%2d/%2d] %-16s UA %-3d  %6.1f s  |  %d -> %d cells (%.0fx smaller)  |  ETA %s",
    i, n_units, sp, ua, dt,
    attr(res, "n_before"), attr(res, "n_after"),
    attr(res, "n_before") / attr(res, "n_after"),
    if (is.na(remaining)) "-" else format(round(remaining / 60, 1), nsmall = 1)
  ))
}

long_metrics <- bind_rows(out_list)

message(sprintf("   Done in %.1f min (%d unit(s) taken from cache).",
                as.numeric(difftime(Sys.time(), t_start, units = "mins")),
                n_cached))


# Map the landscapemetrics metric names onto the column names of the project.
# Forest columns must match Data_Raw_FINAL.csv exactly, because Parts I, II and
# III refer to them by name.
metric_prefix <- c(
  ca      = "CA",
  pland   = "PLAND",
  np      = "NP",
  pd      = "PD",
  ed      = "ED",
  lpi     = "LPI",
  area_mn = "Mean_AREA",
  frac_mn = "Mean_FRAC",
  enn_mn  = "Mean_ENN"
)
metric_suffix <- c(area_mn = "_ha", enn_mn = "_m")   # unit suffixes in use

# Only the metrics actually computed are mapped. With RECOMPUTE_PATCH_MEANS =
# FALSE the three per-patch means are absent, and the agreement check of Section
# 6 simply has fewer columns to compare; the repaired columns are all present.
metric_prefix <- metric_prefix[names(metric_prefix) %in% unique(long_metrics$metric)]

class_labels <- tibble(class = unname(CLASS_CODES), class_name = names(CLASS_CODES))

wide_raster <- long_metrics %>%
  inner_join(class_labels, by = "class") %>%
  mutate(
    column = paste0(metric_prefix[metric], "_", class_name,
                    dplyr::coalesce(metric_suffix[metric], ""))
  ) %>%
  select(SPECIES, UA_ID, UA_Area_ha_recomp, column, value) %>%
  pivot_wider(names_from = column, values_from = value)

message(sprintf("   %d units by %d metric columns.",
                nrow(wide_raster), ncol(wide_raster) - 3))


# ==============================================================================
# 5. PROXIMITY INDEX FROM THE FRAGSTATS PATCH TABLES                   [3/5]
# ==============================================================================
# PROX_MN is the arithmetic mean of the patch-level Proximity Index over the
# patches of the class. For the farming class the patches of cls_3 and cls_4 are
# pooled before the mean is taken, which is the operation that produced the
# stored values for the other 65 units.
#
# This is not the Proximity Index of the merged class in the strict sense: the
# search neighbourhoods were evaluated on the unmerged map. It is, however, the
# quantity that the manuscript reports, and applying it identically to all 67
# units is what keeps the column internally consistent. State this in the
# Methods as the operational definition of PROX_MN for the farming class.
# ==============================================================================
message("\n[3/5] Recovering PROX_MN from the FRAGSTATS patch tables...")

patch_files <- list.files(fragstats_path, pattern = "Patch\\.CSV$",
                          ignore.case = TRUE, full.names = TRUE)
if (length(patch_files) == 0) {
  stop("No *_Patch.CSV found in ", fragstats_path)
}

read_fragstats <- function(f) {
  readr::read_csv(f, col_types = readr::cols(.default = "c"),
                  name_repair = "unique_quiet") %>%
    rename_with(str_squish) %>%
    mutate(across(where(is.character), str_squish))
}

prox_table <- patch_files %>%
  map_dfr(read_fragstats) %>%
  select(any_of(c("SPECIES", "UA", "TYPE", "PROX"))) %>%
  mutate(
    SPECIES = canonical_taxon(SPECIES),
    UA_ID   = as.integer(UA),
    PROX    = as.numeric(PROX),
    class_name = dplyr::case_when(
      TYPE == "cls_1"               ~ "Forest",
      TYPE == "cls_2"               ~ "Herbaceous",
      TYPE %in% c("cls_3", "cls_4") ~ "Agropecuaria",
      TYPE == "cls_5"               ~ "Water",
      TYPE == "cls_6"               ~ "Non_Vegetated",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(SPECIES), !is.na(UA_ID), !is.na(class_name), !is.na(PROX)) %>%
  group_by(SPECIES, UA_ID, class_name) %>%
  summarise(PROX_MN = mean(PROX), .groups = "drop") %>%
  mutate(column = paste0("PROX_MN_", class_name)) %>%
  select(SPECIES, UA_ID, column, PROX_MN) %>%
  pivot_wider(names_from = column, values_from = PROX_MN)

message(sprintf("   PROX_MN recovered for %d units.", nrow(prox_table)))

recomputed <- wide_raster %>%
  left_join(prox_table, by = c("SPECIES", "UA_ID")) %>%
  arrange(SPECIES, UA_ID)

if (nrow(recomputed) != N_UNITS_EXPECTED) {
  stop(sprintf("Recomputed table has %d rows, expected %d.",
               nrow(recomputed), N_UNITS_EXPECTED))
}

write_csv(recomputed, file.path(out_path, "Tabela_Metricas_Recomputada.csv"))
message(sprintf("   Tabela_Metricas_Recomputada.csv exported (%d units).",
                nrow(recomputed)))


# ==============================================================================
# 6. AGREEMENT WITH THE STORED VALUES                                  [4/5]
# ==============================================================================
# Compares every recomputed cell against the corresponding stored cell of
# Data_Raw_FINAL.csv, for the cells where a stored value exists. Cells that are
# empty in the stored file are the target of the repair and are reported
# separately in Section 7.
# ==============================================================================
message("\n[4/5] Comparing recomputed against stored values...")

stored <- read_csv(file.path(out_path, "Data_Raw_FINAL.csv"), show_col_types = FALSE) %>%
  mutate(SPECIES = canonical_taxon(SPECIES),
         UA_ID   = as.integer(UA_ID))

shared_cols <- intersect(
  setdiff(names(recomputed), c("SPECIES", "UA_ID", "UA_Area_ha_recomp")),
  names(stored)
)
message(sprintf("   %d metric columns present in both tables.", length(shared_cols)))

comparison <- shared_cols %>%
  map_dfr(function(cl) {
    j <- stored %>%
      select(SPECIES, UA_ID, stored_value = all_of(cl)) %>%
      inner_join(recomputed %>% select(SPECIES, UA_ID, new_value = all_of(cl)),
                 by = c("SPECIES", "UA_ID")) %>%
      filter(!is.na(stored_value), !is.na(new_value))

    if (nrow(j) < 3) {
      return(tibble(Column = cl, N = nrow(j), Pearson_r = NA_real_,
                    Mean_Bias = NA_real_, Median_Rel_Error = NA_real_,
                    Max_Rel_Error = NA_real_))
    }

    rel <- abs(j$new_value - j$stored_value) /
      ifelse(abs(j$stored_value) < 1e-9, NA_real_, abs(j$stored_value))

    tibble(
      Column           = cl,
      N                = nrow(j),
      Pearson_r        = suppressWarnings(cor(j$stored_value, j$new_value)),
      Mean_Bias        = mean(j$new_value - j$stored_value),
      Median_Rel_Error = median(rel, na.rm = TRUE),
      Max_Rel_Error    = suppressWarnings(max(rel, na.rm = TRUE))
    )
  }) %>%
  mutate(
    Verdict = dplyr::case_when(
      is.na(Median_Rel_Error)     ~ "too few pairs to judge",
      Median_Rel_Error < 0.01     ~ "reproduces the stored pipeline; fill gaps only",
      Median_Rel_Error < 0.05     ~ "small systematic offset; inspect before use",
      TRUE                        ~ "not interchangeable; replace the whole column"
    )
  ) %>%
  arrange(desc(Median_Rel_Error))

write_csv(comparison, file.path(report_path, "TableS_Recompute_Agreement.csv"))

message("\n   Agreement, worst ten metrics first:")
print(as.data.frame(comparison %>% slice_head(n = 10)), row.names = FALSE, digits = 4)

n_bad <- sum(comparison$Median_Rel_Error >= 0.01, na.rm = TRUE)
if (n_bad > 0) {
  message(sprintf(
    "\n   WARNING. %d column(s) do not reproduce the stored values to within 1%%.",
    n_bad))
  message("   For those columns, replace the values of ALL 67 units with the")
  message("   recomputed ones rather than filling only the two gaps, so that a")
  message("   single measurement procedure underlies the column, and record the")
  message("   change in the Methods.")
} else {
  message("\n   All shared columns reproduce the stored values to within 1%.")
  message("   Filling only the empty cells is therefore methodologically safe.")
}


# ==============================================================================
# 7. INVENTORY OF THE CELLS THAT WILL BE REPAIRED                      [5/5]
# ==============================================================================
message("\n[5/5] Listing the cells that Part I will fill...")

repaired <- shared_cols %>%
  map_dfr(function(cl) {
    stored %>%
      select(SPECIES, UA_ID, stored_value = all_of(cl)) %>%
      inner_join(recomputed %>% select(SPECIES, UA_ID, new_value = all_of(cl)),
                 by = c("SPECIES", "UA_ID")) %>%
      filter(is.na(stored_value), !is.na(new_value)) %>%
      transmute(SPECIES, UA_ID, Column = cl, Filled_value = new_value)
  }) %>%
  arrange(SPECIES, UA_ID, Column)

write_csv(repaired, file.path(report_path, "TableS_Repaired_Cells.csv"))

message(sprintf("   %d empty cells will be filled, in %d sampling unit(s):",
                nrow(repaired), n_distinct(paste(repaired$SPECIES, repaired$UA_ID))))
print(as.data.frame(repaired %>% count(SPECIES, UA_ID, name = "cells_filled")),
      row.names = FALSE)

still_empty <- shared_cols %>%
  map_dfr(function(cl) {
    stored %>%
      select(SPECIES, UA_ID, stored_value = all_of(cl)) %>%
      inner_join(recomputed %>% select(SPECIES, UA_ID, new_value = all_of(cl)),
                 by = c("SPECIES", "UA_ID")) %>%
      filter(is.na(stored_value), is.na(new_value)) %>%
      transmute(SPECIES, UA_ID, Column = cl)
  })

if (nrow(still_empty) > 0) {
  message(sprintf(
    "\n   %d cell(s) remain empty after recomputation. These are genuinely",
    nrow(still_empty)))
  message("   undefined, not missing: a per-patch mean over a class that has no")
  message("   patch in the unit, or ENN over a single patch. They are structural")
  message("   zeros only for CA, PLAND, NP, PD, ED and LPI, and Part I recodes")
  message("   them accordingly. They are NOT zeros for AREA_MN, FRAC_MN, ENN_MN")
  message("   or PROX_MN.")
  print(as.data.frame(still_empty %>% count(Column, name = "n_units")),
        row.names = FALSE)
}

message("\n============================================================")
message(sprintf("  Part 0 complete. %d units recomputed.", nrow(recomputed)))
if (!RECOMPUTE_PATCH_MEANS) {
  message("  AREA_MN, FRAC_MN and ENN_MN were not recomputed. They are complete in")
  message("  Data_Raw_FINAL.csv for all 67 units and need no repair. To extend the")
  message("  agreement table to them, set RECOMPUTE_PATCH_MEANS <- TRUE and rerun;")
  message("  the cache keeps the two metric sets separate, so nothing is lost.")
}
message("  Next: run Part I, which reads Tabela_Metricas_Recomputada.csv and")
message("  repairs Data_Raw_FINAL.csv before any analysis.")
message("============================================================\n")

message("\n[Session]")
message("\n[Session] Package versions and system information:")
print(sessionInfo())
writeLines(capture.output(sessionInfo()),
           file.path(report_path, "session_info_part00.txt"))

# ==============================================================================
# END OF PART 0 SCRIPT (v2)
# ==============================================================================
