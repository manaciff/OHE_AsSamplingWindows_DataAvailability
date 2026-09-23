# ==============================================================================
# 01  Data preparation and unconstrained ordination
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
#   Assemble the analytical dataset and describe the landscape inside each
#   occupied habitat envelope. Taxon labels are canonicalized, the metric gaps
#   left by script 00 are repaired, structural zeros are recoded conditionally,
#   and polygon centroids are extracted. Composition is then ordinated by
#   non-metric multidimensional scaling on Bray-Curtis dissimilarities and
#   tested with betadisper before PERMANOVA, because the design is unbalanced.
#
# INPUT
#   Dados/Brutos raw table; Tabela_Metricas_Recomputada.csv; per-species
#   shapefiles
#
# OUTPUT
#   Dados/Processados/Data_Raw_WithCoords.csv, the sole input of scripts 02, 03
#   and 05; Tables 1 to 4; Tables S2 to S5; TableS_NMDS_Summary.csv;
#   TableS_Metric_Repair_Log.csv; TableS_Sample_Sizes_PartI.csv; Figures 1 to 3
#
# GUARDS
#   Stops unless the input has 67 rows, unless all 67 carry coordinates, and
#   unless all 67 are complete over the analysis pool.
#
# FIGURE TITLES
#   Both target journals require the title in the caption rather than inside the
#   image, so every labs(title=) and labs(subtitle=) feeding a manuscript figure
#   was removed and the quantities they carried moved to the captions.
# ==============================================================================


# ==============================================================================
# 1. ENVIRONMENT SETUP
# ==============================================================================
suppressPackageStartupMessages({
  library(here)          # reproducible relative paths
  library(sf)            # spatial vector data (centroid extraction)
  library(tidyverse)     # data wrangling and ggplot2
  library(vegan)         # NMDS, PERMANOVA, betadisper, envfit
  library(ggridges)      # ridgeline density plots
  library(scales)        # axis formatting
  library(corrplot)      # collinearity plots
})

# Reproducibility
GLOBAL_SEED <- 123
N_PERM      <- 9999

# The full sampling design. Every analysis in Parts I, II and III must run on
# this many units; the guards below enforce it.
N_UNITS_EXPECTED <- 67

# Paths (relative to project root)
data_path   <- here("Dados", "Processados")
output_path <- here("Outputs", "Manuscrito", "PartI")
geo_path    <- here("08_Dados_Especies", "Dados_geo_especies")
shp_path    <- file.path(geo_path, "Sp_data_singlepart")

if (!dir.exists(output_path)) dir.create(output_path, recursive = TRUE)


# ==============================================================================
# 1.1 COLORBLIND-SAFE PALETTES (Okabe-Ito, Springer-friendly)
# ==============================================================================
colors_genus <- c(
  "Puma"         = "#D55E00",
  "Alouatta"     = "#009E73",
  "Mazama"       = "#CC79A7",
  "Brachyteles"  = "#E69F00",
  "Tayassu"      = "#0072B2",
  "Myrmecophaga" = "#F0E442",
  "Leopardus"    = "#56B4E9",
  "Bradypus"     = "#999999"
)


# ==============================================================================
# 2. UTILITY FUNCTIONS
# ==============================================================================

save_publication_plot <- function(plot_obj, file_name,
                                  width_mm  = 174,
                                  height_mm = 130,
                                  dpi       = 600) {
  if (is.null(plot_obj)) return(invisible(NULL))
  ggsave(filename    = file_name,
         plot        = plot_obj,
         path        = output_path,
         width       = width_mm,
         height      = height_mm,
         units       = "mm",
         dpi         = dpi,
         device      = "tiff",
         compression = "lzw",
         bg          = "white")
  message(sprintf("   Figure exported: %s (%d x %d mm, %d dpi)",
                  file_name, width_mm, height_mm, dpi))
}

# ------------------------------------------------------------------------------
# canonical_taxon(): one label per taxon, whatever spelling reaches the function.
#
# The project accumulated four kinds of deviation, all of them silent:
#   aguriba        raster ud6_aguriba1.tif, missing the second 'a' of guariba
#   baracnoides    missing the 'h' of arachnoides, in the intermediate table
#   concolor       Puma rasters ud10-ud12 drop the 'p' prefix used by ud1-ud9
#   B_hypoxantus   invalid spelling of Brachyteles hypoxanthus (Kuhl, 1820);
#                  the tabular data used it, the shapefile used the valid one
#
# The function strips every non-letter and lowercases before matching, so the
# form of the separator never matters. It returns NA for an unrecognised name,
# which the callers turn into an error. Silent non-matching is what produced
# n = 65, so an unknown taxon must stop the run, never pass through.
#
# This definition is reproduced verbatim in Part 0. Keep the two in step.
# ------------------------------------------------------------------------------
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
         paste(unique(x[is.na(out) & !is.na(x)]), collapse = ", "),
         "\nAdd them to canonical_taxon() in Part 0 and Part I before rerunning.")
  }
  out
}

tidy_adonis <- function(adonis_obj, group_label) {
  row <- as.data.frame(adonis_obj)[1L, , drop = FALSE]
  data.frame(
    Grouping = group_label,
    Df       = row$Df,
    SumOfSqs = round(row$SumOfSqs, 4),
    R2       = round(row$R2, 4),
    F_value  = round(row$F, 4),
    p_value  = row[["Pr(>F)"]]
  )
}

pairwise.adonis <- function(x, factors,
                            sim.method = "bray",
                            p.adjust.m = "fdr",
                            n_perm     = N_PERM,
                            seed       = GLOBAL_SEED) {
  co <- combn(unique(as.character(factors)), 2)
  pairs <- character(0); F.Model <- numeric(0)
  R2 <- numeric(0); p.value <- numeric(0)

  for (elem in seq_len(ncol(co))) {
    sub_idx <- factors %in% c(co[1, elem], co[2, elem])
    x1 <- vegdist(x[sub_idx, , drop = FALSE], method = sim.method)
    set.seed(seed)
    ad <- adonis2(x1 ~ factors[sub_idx], permutations = n_perm)
    pairs   <- c(pairs, paste(co[1, elem], "vs", co[2, elem]))
    F.Model <- c(F.Model, ad$F[1])
    R2      <- c(R2, ad$R2[1])
    p.value <- c(p.value, ad$`Pr(>F)`[1])
  }

  data.frame(
    pairs      = pairs,
    F.Model    = round(F.Model, 4),
    R2         = round(R2, 4),
    p.value    = round(p.value, 4),
    p.adjusted = round(p.adjust(p.value, method = p.adjust.m), 4),
    sig        = ifelse(p.adjust(p.value, method = p.adjust.m) < 0.05, "*", "")
  )
}


# ==============================================================================
# 3. DATA LOADING & PREPROCESSING                                      [1/9]
# ==============================================================================
# Input: Data_Raw_FINAL.csv — corrected UA_IDs and areas; configuration
# predictors joined from UA_Variables_GLM.xlsx by compound key (SPECIES, UA_ID).
# ==============================================================================
message("\n[1/9] Loading and formatting data...")

data_raw <- read_csv(file.path(data_path, "Data_Raw_FINAL.csv"),
                     show_col_types = FALSE)

# Backwards compatibility for column rename
if ("UA_Area_ha" %in% names(data_raw) &&
    !("AOH_unit_area_ha" %in% names(data_raw))) {
  data_raw <- data_raw %>% rename(AOH_unit_area_ha = UA_Area_ha)
  message("   Column UA_Area_ha renamed to AOH_unit_area_ha.")
}

# Canonical taxon labels, applied before anything is joined or counted.
data_raw$SPECIES <- canonical_taxon(data_raw$SPECIES)
data_raw$UA_ID   <- as.integer(data_raw$UA_ID)

if (nrow(data_raw) != N_UNITS_EXPECTED) {
  stop(sprintf("Data_Raw_FINAL.csv holds %d rows, expected %d.",
               nrow(data_raw), N_UNITS_EXPECTED))
}


# ------------------------------------------------------------------------------
# 3.1 RESPONSE AT FULL PRECISION, WHEN IT IS AVAILABLE
# ------------------------------------------------------------------------------
# ADDED 23 Aug 2026. The area column of Data_Raw_FINAL.csv is stored as a whole
# number: 34 of the 67 values are exact multiples of ten, 28 of a hundred and
# five of a thousand. An area computed by a geographic information system does
# not produce that pattern, so the response was rounded somewhere between ArcGIS
# and the consolidated table.
#
# The consequence is small but it touches the central argument. Part V reports
# the identity log(OHE) = log(NP) + log(100) - log(PD) as closing to within
# 0.029 on the logarithmic scale. Some of that residual is unavoidable, because
# the metrics are counted over 30 m pixels while the response is the area of a
# polygon, and the two do not describe exactly the same surface. The rest is the
# rounding of the response, and that part can be removed.
#
# arcpy_AOH_area_full_precision.py, in the root of the repository, exports the
# area of every sampling unit from the same polygons at full precision. If the
# file it writes is present, it is used here; if it is absent, the run continues
# on the stored values and says so. Nothing is imputed either way.
# ------------------------------------------------------------------------------
precision_file <- file.path(data_path, "AOH_area_full_precision.csv")
if (file.exists(precision_file)) {
  prec <- read_csv(precision_file, show_col_types = FALSE) %>%
    mutate(SPECIES = canonical_taxon(SPECIES),
           UA_ID   = as.integer(UA_ID)) %>%
    dplyr::select(SPECIES, UA_ID, area_planar_ha)

  if (anyDuplicated(prec[, c("SPECIES", "UA_ID")]) > 0)
    stop("AOH_area_full_precision.csv carries duplicated (SPECIES, UA_ID) keys.")

  joined <- data_raw %>%
    dplyr::select(SPECIES, UA_ID, stored = AOH_unit_area_ha) %>%
    left_join(prec, by = c("SPECIES", "UA_ID"))

  missing_units <- sum(is.na(joined$area_planar_ha))
  if (missing_units > 0) {
    print(as.data.frame(joined[is.na(joined$area_planar_ha),
                               c("SPECIES", "UA_ID")]), row.names = FALSE)
    stop(sprintf(paste0("%d unit(s) have no area in AOH_area_full_precision.csv. ",
                        "Rerun arcpy_AOH_area_full_precision.py over all nine ",
                        "shapefiles, or delete the file to run on the stored ",
                        "values."), missing_units))
  }

  rel <- abs(joined$area_planar_ha - joined$stored) / joined$stored
  # A departure of more than five per cent is not rounding; it means the polygons
  # exported are not the ones the stored table describes.
  if (max(rel) > 0.05) {
    worst <- joined[which.max(rel), ]
    stop(sprintf(paste0("Full-precision area departs from the stored value by ",
                        "%.1f%% on %s UA %d (%.1f against %.1f ha). That is too ",
                        "large to be rounding. Check that the shapefiles are the ",
                        "ones the analysis used before replacing the response."),
                 100 * max(rel), worst$SPECIES, worst$UA_ID,
                 worst$area_planar_ha, worst$stored))
  }

  data_raw <- data_raw %>%
    left_join(prec, by = c("SPECIES", "UA_ID")) %>%
    mutate(AOH_unit_area_ha = area_planar_ha) %>%
    dplyr::select(-area_planar_ha)

  message(sprintf(paste0("   Response replaced by the full-precision area: ",
                         "median change %.3f%%, maximum %.3f%%."),
                  100 * median(rel), 100 * max(rel)))
} else {
  message(paste0("   AOH_area_full_precision.csv not found; using the stored ",
                 "response. See arcpy_AOH_area_full_precision.py in the root of ",
                 "the repository."))
}


# ==============================================================================
# 3.5 METRIC REPAIR AND CONDITIONAL STRUCTURAL ZEROS
# ==============================================================================
# Two sampling units, A_guariba UA 6 and B_arachnoides UA 1, entered
# Data_Raw_FINAL.csv with empty cells in ED and LPI of every class, and in NP, PD
# and PROX_MN of the farming class. They were inserted by a reconstruction that
# could only read the exported FRAGSTATS tables, and those tables do not carry
# patch perimeter, nor patch counts for the union of cls_3 and cls_4.
#
# The cells are filled from Tabela_Metricas_Recomputada.csv, produced by
# Script_Part_0_Recompute_Units_v2.R. Nothing is imputed. Each value is measured
# again from the same per-unit raster, with the same reclassification and the
# same 8-neighbour rule that produced the other 65 units, and PROX_MN is taken
# from the same FRAGSTATS patch tables by the same pooled mean. Part 0 exports an
# agreement table quantifying, metric by metric, how closely the recomputation
# reproduces the stored values; that table belongs in the supplementary material.
#
# The distinction that follows is the one the previous pipeline collapsed.
# FRAGSTATS omits the row of a class that is absent from a unit, and the
# consolidation turned that omission into NA. A unit with no herbaceous patch has
# a class area of zero, a percentage of landscape of zero, no patch and no edge;
# it does not have an unknown value. But a unit with 91.6 per cent forest cover
# and an empty ED_Forest cell has a missing value, not a zero. Recoding every NA
# to zero conflated the two and injected two false zeros of forest edge density
# into the variation partitioning and the RDA of Part II.
#
# A cell is therefore zeroed only when the class area of its class is zero. If an
# NA survives in a class that is present, the run stops and names the cells,
# because that means Part 0 was not run, or did not cover the unit.
# ==============================================================================
message("\n[1.5] Repairing metric gaps and recoding structural zeros...")

recomp_file <- file.path(data_path, "Tabela_Metricas_Recomputada.csv")
if (!file.exists(recomp_file)) {
  stop("Tabela_Metricas_Recomputada.csv not found in: ", data_path, "\n",
       "Run Script_Part_0_Recompute_Units_v2.R first. Without it, two sampling ",
       "units carry empty ED, LPI and farming-class cells, and every analysis ",
       "downstream silently drops them.")
}

recomp <- read_csv(recomp_file, show_col_types = FALSE) %>%
  mutate(SPECIES = canonical_taxon(SPECIES),
         UA_ID   = as.integer(UA_ID))

repair_cols <- intersect(
  setdiff(names(recomp), c("SPECIES", "UA_ID", "UA_Area_ha_recomp")),
  names(data_raw)
)

n_before <- nrow(data_raw)
data_raw <- data_raw %>%
  left_join(
    recomp %>%
      select(SPECIES, UA_ID, all_of(repair_cols)) %>%
      rename_with(~paste0(.x, "__recomp"), all_of(repair_cols)) %>%
      # Explicit marker of a successful match. Testing the first repair column
      # instead would be unreliable, because a column such as LPI_Water is
      # legitimately NA wherever the class is absent, and Part 0 may export a
      # reduced metric set when RECOMPUTE_PATCH_MEANS is FALSE, which changes
      # which column comes first.
      mutate(matched__recomp = TRUE),
    by = c("SPECIES", "UA_ID")
  )
stopifnot(nrow(data_raw) == n_before)

n_unmatched <- sum(is.na(data_raw$matched__recomp))
if (n_unmatched > 0) {
  print(as.data.frame(
    data_raw[is.na(data_raw$matched__recomp), c("SPECIES", "UA_ID")]),
    row.names = FALSE)
  stop(sprintf(paste0("%d unit(s) of Data_Raw_FINAL.csv have no counterpart in ",
                      "the recomputed table. The join key is (SPECIES, UA_ID) ",
                      "after canonical_taxon(); check Part 0."), n_unmatched))
}
data_raw$matched__recomp <- NULL

message(sprintf("   Recomputed table joined: %d column(s) available for repair.",
                length(repair_cols)))

repair_log <- list()
n_filled   <- 0L
for (cl in repair_cols) {
  src <- paste0(cl, "__recomp")
  idx <- is.na(data_raw[[cl]]) & !is.na(data_raw[[src]])
  if (any(idx)) {
    repair_log[[length(repair_log) + 1L]] <- data.frame(
      SPECIES = data_raw$SPECIES[idx],
      UA_ID   = data_raw$UA_ID[idx],
      Column  = cl,
      Value   = data_raw[[src]][idx]
    )
    data_raw[[cl]][idx] <- data_raw[[src]][idx]
    n_filled <- n_filled + sum(idx)
  }
}
data_raw <- data_raw %>% select(-ends_with("__recomp"))

if (n_filled > 0) {
  repair_df <- bind_rows(repair_log)
  write.csv(repair_df,
            file.path(output_path, "TableS_Metric_Repair_Log.csv"),
            row.names = FALSE)
  message(sprintf("   %d cell(s) filled from the recomputed table, in %d unit(s):",
                  n_filled, n_distinct(paste(repair_df$SPECIES, repair_df$UA_ID))))
  print(as.data.frame(repair_df %>% count(SPECIES, UA_ID, name = "cells_filled")),
        row.names = FALSE)
} else {
  message("   No gap to fill: Data_Raw_FINAL.csv is already complete.")
}

# --- conditional structural zeros -------------------------------------------
CLASSES_LULC <- c("Forest", "Herbaceous", "Agropecuaria", "Water", "Non_Vegetated")
ZERO_METRICS <- c("CA", "PLAND", "NP", "PD", "ED", "LPI")

n_zeroed  <- 0L
bad_cells <- list()

for (cls in CLASSES_LULC) {
  ca_col <- paste0("CA_", cls)
  if (!ca_col %in% names(data_raw)) next

  ca_val    <- data_raw[[ca_col]]
  is_absent <- is.na(ca_val) | ca_val <= 0

  for (mt in ZERO_METRICS) {
    col <- paste0(mt, "_", cls)
    if (!col %in% names(data_raw)) next

    to_zero <- is_absent & is.na(data_raw[[col]])
    if (any(to_zero)) {
      data_raw[[col]][to_zero] <- 0
      n_zeroed <- n_zeroed + sum(to_zero)
    }

    still_na <- !is_absent & is.na(data_raw[[col]])
    if (any(still_na)) {
      bad_cells[[length(bad_cells) + 1L]] <- data.frame(
        SPECIES = data_raw$SPECIES[still_na],
        UA_ID   = data_raw$UA_ID[still_na],
        Class   = cls,
        Column  = col,
        CA      = ca_val[still_na]
      )
    }
  }
}

message(sprintf("   %d structural zeros recoded (class genuinely absent).", n_zeroed))

# The stop is scoped to the three focal classes. Water and Non-Vegetated were
# excluded at the KDE stage and enter no analysis, so an unresolved cell in them
# is worth reporting but is not a reason to halt the pipeline.
FOCAL_CLASSES <- c("Forest", "Herbaceous", "Agropecuaria")

if (length(bad_cells) > 0) {
  bad_df   <- bind_rows(bad_cells)
  bad_focal <- bad_df %>% filter(Class %in% FOCAL_CLASSES)
  bad_other <- bad_df %>% filter(!Class %in% FOCAL_CLASSES)

  if (nrow(bad_other) > 0) {
    message(sprintf("   %d cell(s) unresolved in non-focal classes (not used in any analysis):",
                    nrow(bad_other)))
    print(as.data.frame(bad_other %>% count(Column, name = "n_units")),
          row.names = FALSE)
  }

  if (nrow(bad_focal) > 0) {
    print(as.data.frame(bad_focal), row.names = FALSE)
    stop("The cells above are NA in a focal class that is present in the unit. ",
         "They are missing values, not structural zeros, and must not be ",
         "recoded to zero. Run Script_Part_0_Recompute_Units_v2.R.")
  }
}

# Per-patch means are left as NA where they are genuinely undefined: a mean over
# zero patches, or a nearest-neighbour distance for a single patch. These are not
# zeros and are not filled. They are reported so that the count is on record.
undefined_cols <- grep("^(PROX_MN|Mean_ENN|Mean_FRAC|Mean_AREA)_",
                       names(data_raw), value = TRUE)
n_undef <- sum(is.na(data_raw[, undefined_cols]))
if (n_undef > 0) {
  message(sprintf("   %d per-patch mean(s) remain undefined (no patch of the class).",
                  n_undef))
  undef_report <- data_raw %>%
    select(SPECIES, UA_ID, all_of(undefined_cols)) %>%
    pivot_longer(-c(SPECIES, UA_ID), names_to = "Column", values_to = "Value") %>%
    filter(is.na(Value)) %>%
    count(Column, name = "n_units")
  print(as.data.frame(undef_report), row.names = FALSE)
}


# ==============================================================================
# 3.6 GUARD ON THE ANALYTICAL SAMPLE SIZE
# ==============================================================================
# Every analysis of Parts I, II and III draws its predictors from the pool below.
# If a single unit is incomplete over that pool, complete.cases() in Part II and
# drop_na() in Part III will remove it, and the manuscript will report a sample
# size smaller than the design. That is exactly how n = 63 and n = 65 reached the
# supplementary tables. The run stops here instead, naming the offending cells.
# ==============================================================================
ANALYSIS_POOL <- c(
  "PLAND_Forest", "PLAND_Agropecuaria", "PLAND_Herbaceous",
  "PD_Forest", "ED_Forest", "Mean_FRAC_Forest", "Mean_ENN_Forest_m",
  "PROX_MN_Forest", "Mean_AREA_Forest_ha",
  "PD_Agropecuaria", "ED_Agropecuaria", "PROX_MN_Agropecuaria",
  "Mean_Elevation_m", "Mean_Pop_Density"
)
pool_present <- intersect(ANALYSIS_POOL, names(data_raw))
n_complete   <- sum(complete.cases(data_raw[, pool_present]))

message(sprintf("   Units complete over the %d-predictor analysis pool: %d of %d.",
                length(pool_present), n_complete, nrow(data_raw)))

if (n_complete != N_UNITS_EXPECTED) {
  offenders <- data_raw %>%
    filter(!complete.cases(data_raw[, pool_present])) %>%
    select(SPECIES, UA_ID, all_of(pool_present)) %>%
    pivot_longer(-c(SPECIES, UA_ID), names_to = "Column", values_to = "Value") %>%
    filter(is.na(Value))
  print(as.data.frame(offenders), row.names = FALSE)
  stop(sprintf(paste0("Only %d of %d units are complete over the analysis pool. ",
                      "Parts II and III would drop the rest silently. Resolve ",
                      "the cells listed above before proceeding."),
               n_complete, N_UNITS_EXPECTED))
}


# Factor conversions
data_raw <- data_raw %>%
  mutate(across(c(ORDER, GENUS, SPECIES, DIET, LOCOMOTION), as.factor),
         UA_ID = as.integer(UA_ID))

message(sprintf("   Data loaded: %d observations across %d species and %d genera.",
                nrow(data_raw),
                n_distinct(data_raw$SPECIES),
                n_distinct(data_raw$GENUS)))

# Integrity check: (SPECIES, UA_ID) must be unique before spatial join
dup_check <- data_raw %>%
  group_by(SPECIES, UA_ID) %>%
  summarise(n = n(), .groups = "drop") %>%
  filter(n > 1)
if (nrow(dup_check) > 0) {
  warning("Duplicated (SPECIES, UA_ID) pairs detected. Spatial join will be ambiguous.")
  print(dup_check)
} else {
  message("   Integrity OK: every (SPECIES, UA_ID) pair is unique.")
}


# ==============================================================================
# 4. CENTROID EXTRACTION FROM PER-SPECIES SHAPEFILES               (NEW [2/9])
# ==============================================================================
# Reads every .shp file from shp_path and extracts polygon centroids.
# All layers are re-projected to the CRS of the first file for consistency.
#
# Field mapping:
#   Species identifier  -> SP_ID field (preferred) or SPECIES field (fallback).
#                          If neither is present the species name is inferred
#                          from the file stem (CRS suffix stripped).
#   Sampling-unit ID    -> ORIG_FID field (mandatory; execution stops with an
#                          informative error if the field is absent).
#
# Join key: (SPECIES, UA_ID) in data_raw matched against (sp_id, ORIG_FID)
# in the shapefiles — guarantees a one-to-one match.
# ==============================================================================
message("\n[2/9] Extracting polygon centroids from per-species shapefiles...")

shp_files <- list.files(shp_path, pattern = "\\.shp$", full.names = TRUE)
if (length(shp_files) == 0) {
  stop(sprintf(
    "No .shp files found in: %s\nCheck shp_path in Section 1.", shp_path
  ))
}

ref_crs <- st_crs(st_read(shp_files[1], quiet = TRUE))

extract_centroids <- function(file) {
  layer <- st_read(file, quiet = TRUE)
  if (st_crs(layer) != ref_crs) layer <- st_transform(layer, ref_crs)

  nm <- names(layer)

  # Detect species identifier: SP_ID preferred, SPECIES as fallback
  sp_field   <- nm[grep("^(sp_id|species)$", nm, ignore.case = TRUE)][1]

  # Detect sampling-unit identifier: ORIG_FID required
  orig_field <- nm[grep("orig_?fid", nm, ignore.case = TRUE)][1]
  if (is.na(orig_field)) {
    stop(sprintf(
      "ORIG_FID field not found in: %s\nAvailable fields: %s",
      basename(file), paste(nm, collapse = ", ")
    ))
  }

  # Build species vector, then canonicalise. The shapefile of Brachyteles
  # hypoxanthus carries SP_ID = "B_hypoxanthus" while the tabular data carried
  # "B_hypoxantus"; without this step the compound key never matches and the two
  # units of the species reach Part II with no coordinates.
  sp_id <- if (!is.na(sp_field)) {
    as.character(layer[[sp_field]])
  } else {
    rep(
      sub("_(utm|albers|wgs).*$", "",
          tools::file_path_sans_ext(basename(file)),
          ignore.case = TRUE),
      nrow(layer)
    )
  }
  sp_id <- canonical_taxon(sp_id)

  # Invalid rings would make st_centroid() return an empty geometry and drop the
  # unit without a message, so the geometry is repaired first. For a polygon with
  # a concave outline the centroid can also fall outside the polygon; when that
  # happens the representative point on the surface is used instead, which keeps
  # the coordinate inside the sampling unit it describes.
  geom <- st_make_valid(st_geometry(layer))
  cent <- suppressWarnings(st_centroid(geom))

  outside <- !suppressMessages(
    as.logical(sf::st_within(cent, geom, sparse = FALSE)[cbind(seq_along(cent),
                                                               seq_along(cent))])
  )
  outside[is.na(outside)] <- TRUE
  if (any(outside)) {
    cent[outside] <- suppressWarnings(st_point_on_surface(geom[outside]))
    message(sprintf("   %s: %d centroid(s) fell outside the polygon; replaced by st_point_on_surface().",
                    basename(file), sum(outside)))
  }

  xy <- st_coordinates(cent)

  data.frame(
    SPECIES = sp_id,
    UA_ID   = as.integer(layer[[orig_field]]),
    coord_x = xy[, 1],
    coord_y = xy[, 2],
    stringsAsFactors = FALSE
  )
}

ua_coords <- bind_rows(lapply(shp_files, extract_centroids)) %>%
  mutate(SPECIES = canonical_taxon(SPECIES),
         UA_ID   = as.integer(UA_ID))

message(sprintf("   Centroid table built: %d records from %d shapefiles.",
                nrow(ua_coords), length(shp_files)))

# A duplicated key would make the join multiply rows instead of matching them.
if (anyDuplicated(ua_coords[, c("SPECIES", "UA_ID")]) > 0) {
  print(ua_coords[duplicated(ua_coords[, c("SPECIES", "UA_ID")]),
                  c("SPECIES", "UA_ID")])
  stop("Duplicated (SPECIES, UA_ID) in the centroid table.")
}

# Join centroids to data_raw by compound key (SPECIES, UA_ID)
data_raw <- data_raw %>%
  mutate(SPECIES_char = as.character(SPECIES)) %>%
  left_join(ua_coords,
            by = c("SPECIES_char" = "SPECIES", "UA_ID" = "UA_ID")) %>%
  select(-SPECIES_char)

n_matched <- sum(!is.na(data_raw$coord_x))
message(sprintf("   Centroids matched: %d of %d (%.1f%%).",
                n_matched, nrow(data_raw),
                100 * n_matched / nrow(data_raw)))

# A missing centroid is not a tolerable loss. It removes the unit from every
# spatial analysis of Part II, and the previous version only warned, so the loss
# surfaced later as a reduced sample size in a supplementary table.
if (n_matched < nrow(data_raw)) {
  unmatched <- data_raw %>%
    filter(is.na(coord_x)) %>%
    distinct(SPECIES, UA_ID)
  print(as.data.frame(unmatched), row.names = FALSE)
  stop(sprintf(paste0("%d sampling unit(s) have no centroid. Check that the ",
                      "shapefile of the taxon exists in %s and that its SP_ID ",
                      "and ORIG_FID match (SPECIES, UA_ID) after ",
                      "canonical_taxon()."),
               nrow(unmatched), shp_path))
}

stopifnot(n_matched == N_UNITS_EXPECTED)
message(sprintf("   All %d sampling units carry centroid coordinates.", n_matched))


# ==============================================================================
# 5. MULTIVARIATE COMMUNITY ANALYSES                                   [3/9]
# ==============================================================================
# Composition matrix: three focal LULC classes only (Forest, Farming,
# Shrub/Herb). Non-Vegetated and Water are excluded, consistent with their
# removal at the KDE stage.
# LPI_Forest is excluded from the envfit configuration set because it is
# biologically redundant with PLAND_Forest and causes VIF > 40 in Part II.
# ==============================================================================
message("\n[3/9] Multivariate analyses (NMDS, betadisper, PERMANOVA, envfit)...")

# --- 3 focal LULC classes; English labels ---
comp_vars <- data_raw %>%
  select(CA_Forest, CA_Herbaceous, CA_Agropecuaria) %>%
  mutate(across(everything(), ~replace_na(.x, 0)))
names(comp_vars) <- c("Forest", "Shrub/Herb", "Farming")

# Spearman correlation among the 3 focal classes
cor_matrix_comp <- cor(comp_vars, use = "complete.obs", method = "spearman")
tiff(filename    = file.path(output_path, "FigS2_Composition_Correlation.tiff"),
     width       = 174, height = 174, units = "mm",
     res         = 600, compression = "lzw")
corrplot(cor_matrix_comp, method = "number", type = "upper",
         tl.col = "black", tl.srt = 45,
         number.cex = 0.85,
         mar = c(0, 0, 2, 0),
         title = "Spearman correlation among LULC classes (3 focal classes)")
invisible(dev.off())

# NMDS on Bray-Curtis distances of sqrt-transformed Class Area
LULC_sqrt <- sqrt(comp_vars)
LULC_dist <- vegdist(LULC_sqrt, method = "bray", na.rm = TRUE)

set.seed(GLOBAL_SEED)
nmds_result <- metaMDS(LULC_dist, k = 2, trymax = 200, trace = FALSE)
message(sprintf("   NMDS stress (k = 2): %.3f", nmds_result$stress))

# The stress is quoted in the Results and in the caption of the ordination
# figure, so it is written to disk rather than left in the console log.
write.csv(
  data.frame(
    Metric = c("NMDS stress (k = 2)", "Dimensions", "N sampling units",
               "Distance", "Transformation", "Permutations"),
    Value  = c(round(nmds_result$stress, 4), 2, nrow(LULC_sqrt),
               "Bray-Curtis", "square root of Class Area", N_PERM)
  ),
  file.path(output_path, "TableS_NMDS_Summary.csv"), row.names = FALSE
)
message("   TableS_NMDS_Summary.csv exported (stress for the Results and captions).")

# envfit: configuration metrics (LPI_Forest excluded — see revision note iv)
config_metrics <- data_raw %>%
  select(Mean_FRAC_Forest, Mean_ENN_Forest_m, PROX_MN_Forest,
         Mean_AREA_Forest_ha, PLAND_Forest, PD_Forest, ED_Forest)

# envfit(na.rm = TRUE) discards any unit with an NA in the predictor matrix and
# does not report how many. Before the repair of Section 3.5, ED_Forest was NA in
# two units and envfit therefore ran on 65 while the ordination itself ran on 67.
# The count is now asserted, so the two can never diverge again unnoticed.
n_envfit <- sum(complete.cases(config_metrics))
message(sprintf("   envfit fitted on %d of %d units.", n_envfit, nrow(config_metrics)))
if (n_envfit != N_UNITS_EXPECTED) {
  print(as.data.frame(
    data_raw[!complete.cases(config_metrics), c("SPECIES", "UA_ID")]),
    row.names = FALSE)
  stop(sprintf("envfit would run on %d units instead of %d.",
               n_envfit, N_UNITS_EXPECTED))
}

set.seed(GLOBAL_SEED)
envfit_result <- envfit(nmds_result, config_metrics,
                        permutations = N_PERM, na.rm = TRUE)

vetores_df <- as.data.frame(scores(envfit_result, display = "vectors"))
vetores_df$Metric  <- rownames(vetores_df)
vetores_df$R2      <- envfit_result$vectors$r
vetores_df$p_value <- envfit_result$vectors$pvals
write.csv(vetores_df,
          file.path(output_path, "TableS2_Envfit_Configuration.csv"),
          row.names = FALSE)

vetores_sig <- subset(vetores_df, p_value < 0.05)
vetores_sig$NMDS1  <- vetores_sig$NMDS1 * 1.2
vetores_sig$NMDS2  <- vetores_sig$NMDS2 * 1.2
vetores_sig$Metric <- gsub("_Forest", "", vetores_sig$Metric)
message(sprintf("   envfit: %d significant predictors (p < 0.05).",
                nrow(vetores_sig)))

# betadisper -- homogeneity of multivariate dispersion.
# Significance is assessed with permutest() (permutation test) rather than
# anova(): the group distances to the centroid are not normally distributed, so
# the parametric F-test assumed by anova() is not appropriate. permutest()
# evaluates the same observed F statistic against its permutation null
# distribution (Anderson & Walsh 2013). The observed F values are identical to
# those anova() would return; only the p-values and their estimation change.
# A fixed seed is set before each call because the permutation step is stochastic.
disp_genus   <- betadisper(LULC_dist, data_raw$GENUS)
disp_species <- betadisper(LULC_dist, data_raw$SPECIES)
disp_diet    <- betadisper(LULC_dist, data_raw$DIET)
disp_loc     <- betadisper(LULC_dist, data_raw$LOCOMOTION)

set.seed(GLOBAL_SEED); permd_genus   <- permutest(disp_genus,   permutations = N_PERM)
set.seed(GLOBAL_SEED); permd_species <- permutest(disp_species, permutations = N_PERM)
set.seed(GLOBAL_SEED); permd_diet    <- permutest(disp_diet,    permutations = N_PERM)
set.seed(GLOBAL_SEED); permd_loc     <- permutest(disp_loc,     permutations = N_PERM)

betadisper_results <- bind_rows(
  as.data.frame(permd_genus$tab)   %>% mutate(Group = "Genus",      Source = rownames(.)),
  as.data.frame(permd_species$tab) %>% mutate(Group = "Species",    Source = rownames(.)),
  as.data.frame(permd_diet$tab)    %>% mutate(Group = "Diet",       Source = rownames(.)),
  as.data.frame(permd_loc$tab)     %>% mutate(Group = "Locomotion", Source = rownames(.))
) %>% select(Group, Source, everything())

message(sprintf(
  "   betadisper (permutest, %d perm): Genus p = %.4f | Species p = %.4f | Diet p = %.4f | Locomotion p = %.4f",
  N_PERM,
  permd_genus$tab[1, "Pr(>F)"],   permd_species$tab[1, "Pr(>F)"],
  permd_diet$tab[1, "Pr(>F)"],    permd_loc$tab[1, "Pr(>F)"]))

write.csv(betadisper_results,
          file.path(output_path, "TableS3_Betadisper.csv"),
          row.names = FALSE)

# PERMANOVA
set.seed(GLOBAL_SEED)
adonis_genus   <- adonis2(LULC_dist ~ GENUS,      data = data_raw, permutations = N_PERM)
set.seed(GLOBAL_SEED)
adonis_species <- adonis2(LULC_dist ~ SPECIES,    data = data_raw, permutations = N_PERM)
set.seed(GLOBAL_SEED)
adonis_diet    <- adonis2(LULC_dist ~ DIET,       data = data_raw, permutations = N_PERM)
set.seed(GLOBAL_SEED)
adonis_loc     <- adonis2(LULC_dist ~ LOCOMOTION, data = data_raw, permutations = N_PERM)

permanova_global <- bind_rows(
  tidy_adonis(adonis_genus,   "Genus"),
  tidy_adonis(adonis_species, "Species"),
  tidy_adonis(adonis_diet,    "Diet"),
  tidy_adonis(adonis_loc,     "Locomotion")
)
write.csv(permanova_global,
          file.path(output_path, "TableS4_PERMANOVA_Global.csv"),
          row.names = FALSE)

write.csv(pairwise.adonis(LULC_sqrt, data_raw$GENUS),
          file.path(output_path, "TableS5_PERMANOVA_Pairwise_Genus.csv"),
          row.names = FALSE)
write.csv(pairwise.adonis(LULC_sqrt, data_raw$SPECIES),
          file.path(output_path, "TableS5_PERMANOVA_Pairwise_Species.csv"),
          row.names = FALSE)
write.csv(pairwise.adonis(LULC_sqrt, data_raw$DIET),
          file.path(output_path, "TableS5_PERMANOVA_Pairwise_Diet.csv"),
          row.names = FALSE)
write.csv(pairwise.adonis(LULC_sqrt, data_raw$LOCOMOTION),
          file.path(output_path, "TableS5_PERMANOVA_Pairwise_Locomotion.csv"),
          row.names = FALSE)
message("   Pairwise PERMANOVA tables exported.")


# ==============================================================================
# 6. DESCRIPTIVE TABLES (Composition + Configuration by species)      [4/9]
# ==============================================================================
message("\n[4/9] Building descriptive tables by species...")

# Shannon diversity from 3 focal classes (Water and Non-Vegetated excluded)
comp_vars_pland <- data_raw %>%
  select(PLAND_Forest, PLAND_Agropecuaria, PLAND_Herbaceous) %>%
  mutate(across(everything(), ~replace_na(.x, 0)))
data_raw$Shannon_Index <- diversity(comp_vars_pland, index = "shannon")

data_raw$MPS_Forest_ha <- data_raw$CA_Forest / data_raw$NP_Forest
data_raw$MPS_Forest_ha[is.infinite(data_raw$MPS_Forest_ha)] <- NA

# Table 1 — Composition (3 focal LULC classes; English column names)
table_composition <- data_raw %>%
  group_by(SPECIES) %>%
  summarise(
    N_SUs                  = n(),
    PLAND_Forest_Mean      = mean(PLAND_Forest,         na.rm = TRUE),
    PLAND_Forest_SD        = sd(PLAND_Forest,           na.rm = TRUE),
    PLAND_Forest_Median    = median(PLAND_Forest,       na.rm = TRUE),
    PLAND_Forest_Min       = min(PLAND_Forest,          na.rm = TRUE),
    PLAND_Forest_Max       = max(PLAND_Forest,          na.rm = TRUE),
    PLAND_Farming_Mean     = mean(PLAND_Agropecuaria,   na.rm = TRUE),
    PLAND_Farming_SD       = sd(PLAND_Agropecuaria,     na.rm = TRUE),
    PLAND_Farming_Median   = median(PLAND_Agropecuaria, na.rm = TRUE),
    PLAND_ShrubHerb_Mean   = mean(PLAND_Herbaceous,     na.rm = TRUE),
    PLAND_ShrubHerb_SD     = sd(PLAND_Herbaceous,       na.rm = TRUE),
    PLAND_ShrubHerb_Median = median(PLAND_Herbaceous,   na.rm = TRUE),
    Shannon_Mean           = mean(Shannon_Index,        na.rm = TRUE),
    Shannon_SD             = sd(Shannon_Index,          na.rm = TRUE),
    Shannon_Median         = median(Shannon_Index,      na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(PLAND_Forest_Mean))

write.csv(table_composition,
          file.path(output_path, "Table1_Composition_By_Species.csv"),
          row.names = FALSE)

# Table 2 — Configuration (LPI_Forest removed; see revision note iv)
table_configuration <- data_raw %>%
  group_by(SPECIES) %>%
  summarise(
    N_SUs                  = n(),
    NP_Forest_Mean         = mean(NP_Forest,             na.rm = TRUE),
    NP_Forest_Median       = median(NP_Forest,           na.rm = TRUE),
    PD_Forest_Mean         = mean(PD_Forest,             na.rm = TRUE),
    PD_Forest_SD           = sd(PD_Forest,               na.rm = TRUE),
    PD_Forest_Median       = median(PD_Forest,           na.rm = TRUE),
    AREA_MN_Forest_Mean    = mean(Mean_AREA_Forest_ha,   na.rm = TRUE),
    AREA_MN_Forest_SD      = sd(Mean_AREA_Forest_ha,     na.rm = TRUE),
    AREA_MN_Forest_Median  = median(Mean_AREA_Forest_ha, na.rm = TRUE),
    MPS_Forest_Mean        = mean(MPS_Forest_ha,         na.rm = TRUE),
    MPS_Forest_Median      = median(MPS_Forest_ha,       na.rm = TRUE),
    ED_Forest_Mean         = mean(ED_Forest,             na.rm = TRUE),
    ED_Forest_Median       = median(ED_Forest,           na.rm = TRUE),
    FRAC_Forest_Mean       = mean(Mean_FRAC_Forest,      na.rm = TRUE),
    FRAC_Forest_Median     = median(Mean_FRAC_Forest,    na.rm = TRUE),
    ENN_Forest_Mean        = mean(Mean_ENN_Forest_m,     na.rm = TRUE),
    ENN_Forest_Median      = median(Mean_ENN_Forest_m,   na.rm = TRUE),
    PROX_Forest_Mean       = mean(PROX_MN_Forest,        na.rm = TRUE),
    PROX_Forest_Median     = median(PROX_MN_Forest,      na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(AREA_MN_Forest_Mean))

write.csv(table_configuration,
          file.path(output_path, "Table2_Configuration_By_Species.csv"),
          row.names = FALSE)
message("   Composition and Configuration tables exported.")


# ==============================================================================
# 7. FIGURE 1 — NMDS WITH ENVFIT VECTORS (genus level, convex hulls) [5/9]
# ==============================================================================
message("\n[5/9] Building Figure 1 (NMDS + envfit)...")

nmds_scores <- as.data.frame(scores(nmds_result, display = "sites"))
nmds_scores$Genus <- data_raw$GENUS

hulls_genus <- nmds_scores %>%
  group_by(Genus) %>%
  slice(chull(NMDS1, NMDS2)) %>%
  ungroup()

genera_present      <- unique(data_raw$GENUS)
colors_genus_subset <- colors_genus[names(colors_genus) %in% genera_present]
shapes_genus        <- setNames(seq_along(genera_present), genera_present)

plot_nmds_genus <- ggplot() +
  geom_polygon(data  = hulls_genus,
               aes(x = NMDS1, y = NMDS2,
                   fill = Genus, group = Genus),
               alpha = 0.30, colour = "grey30") +
  geom_point(data  = nmds_scores,
             aes(x = NMDS1, y = NMDS2,
                 colour = Genus, shape = Genus),
             size = 2.6, alpha = 0.9) +
  geom_segment(data = vetores_sig,
               aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
               arrow = arrow(length = unit(0.22, "cm")),
               colour = "grey15", linewidth = 0.9, alpha = 0.8) +
  geom_text(data = vetores_sig,
            aes(x = NMDS1 * 1.10, y = NMDS2 * 1.10, label = Metric),
            colour = "black", size = 3.2, fontface = "bold") +
  scale_fill_manual(values  = colors_genus_subset, na.value = "grey90") +
  scale_colour_manual(values = colors_genus_subset, na.value = "grey50") +
  scale_shape_manual(values  = shapes_genus) +
  labs(x = "NMDS 1", y = "NMDS 2") +
  theme_classic(base_size = 10) +
  theme(plot.title    = element_text(face = "bold", size = 11, hjust = 0),
        legend.title  = element_text(face = "bold", size = 9),
        legend.text   = element_text(size = 8),
        legend.position = "right")

save_publication_plot(plot_nmds_genus,
                      file_name = "Figure1_NMDS_Composition_Genus.tiff",
                      width_mm  = 174, height_mm = 120)


# ==============================================================================
# 8. RIDGELINE PLOTS — Landscape characterization within AOH           [6/9]
# ==============================================================================
message("\n[6/9] Building Figure 2 (composition ridgeline) and Figure 3 (configuration ridgeline)...")

# --- Figure 2: LULC composition — 3 focal classes with English labels ---
df_long_comp <- data_raw %>%
  select(GENUS,
         `Forest`     = PLAND_Forest,
         `Farming`    = PLAND_Agropecuaria,
         `Shrub/Herb` = PLAND_Herbaceous) %>%
  pivot_longer(cols      = -GENUS,
               names_to  = "Metric",
               values_to = "Value") %>%
  filter(!is.na(Value)) %>%
  mutate(Metric = factor(Metric,
                         levels = c("Forest", "Farming", "Shrub/Herb")))

median_comp <- df_long_comp %>%
  group_by(Metric, GENUS) %>%
  summarise(Median = median(Value, na.rm = TRUE), .groups = "drop")

genus_order_comp <- df_long_comp %>%
  filter(Metric == "Forest") %>%
  group_by(GENUS) %>%
  summarise(med = median(Value, na.rm = TRUE), .groups = "drop") %>%
  arrange(med) %>%
  pull(GENUS) %>% as.character()

df_long_comp$GENUS <- factor(df_long_comp$GENUS, levels = genus_order_comp)
median_comp$GENUS  <- factor(median_comp$GENUS,  levels = genus_order_comp)

plot_ridge_comp <- ggplot(df_long_comp,
                          aes(x = Value, y = GENUS, fill = GENUS)) +
  geom_density_ridges(scale            = 1.15,
                      rel_min_height   = 0.01,
                      alpha            = 0.75,
                      colour           = "grey25",
                      linewidth        = 0.35,
                      quantile_lines   = TRUE,
                      quantiles        = 0.5,
                      vline_colour     = "grey10",
                      vline_width      = 0.8) +
  scale_fill_manual(values = colors_genus[genus_order_comp]) +
  scale_x_continuous(labels = label_number(suffix = "%", accuracy = 1),
                     expand = expansion(mult = c(0.02, 0.05))) +
  facet_wrap(~ Metric, scales = "free_x", nrow = 2) +
  labs(x = "Percentage of landscape (PLAND)", y = NULL) +
  theme_ridges(grid = TRUE, center_axis_labels = TRUE) +
  theme(legend.position    = "none",
        strip.text         = element_text(face = "bold", size = 10,
                                          colour = "white"),
        strip.background   = element_rect(fill = "grey25", colour = NA),
        axis.text.y        = element_text(face = "italic", size = 9),
        axis.text.x        = element_text(size = 8),
        axis.title.x       = element_text(face = "bold", size = 10,
                                          margin = margin(t = 6)),
        panel.spacing      = unit(1, "lines"),
        plot.margin        = margin(6, 8, 6, 6))

save_publication_plot(plot_ridge_comp,
                      file_name = "Figure2_Ridgeline_Composition.tiff",
                      width_mm  = 174, height_mm = 150)

# --- Figure 3: Forest configuration metrics ---
df_long_conf <- data_raw %>%
  transmute(
    GENUS,
    `PD (n / 100 ha)` = PD_Forest,
    `AREA_MN (ha)`    = Mean_AREA_Forest_ha,
    `ED (m / ha)`     = ED_Forest,
    `FRAC_MN`         = Mean_FRAC_Forest,
    `ENN_MN (m)`      = Mean_ENN_Forest_m,
    `PROX_MN`         = PROX_MN_Forest
  ) %>%
  pivot_longer(cols      = -GENUS,
               names_to  = "Metric",
               values_to = "Value") %>%
  filter(!is.na(Value), is.finite(Value)) %>%
  mutate(Metric = factor(Metric,
                         levels = c("PD (n / 100 ha)", "AREA_MN (ha)",
                                    "ED (m / ha)",     "FRAC_MN",
                                    "ENN_MN (m)",      "PROX_MN")))

df_long_conf$GENUS <- factor(df_long_conf$GENUS, levels = genus_order_comp)

plot_ridge_conf <- ggplot(df_long_conf,
                          aes(x = Value, y = GENUS, fill = GENUS)) +
  geom_density_ridges(scale            = 1.10,
                      rel_min_height   = 0.01,
                      alpha            = 0.75,
                      colour           = "grey25",
                      linewidth        = 0.35,
                      quantile_lines   = TRUE,
                      quantiles        = 0.5,
                      vline_colour     = "grey10",
                      vline_width      = 0.8) +
  scale_fill_manual(values = colors_genus[genus_order_comp]) +
  scale_x_continuous(labels = label_number(accuracy = 0.1, big.mark = ","),
                     expand = expansion(mult = c(0.02, 0.05))) +
  facet_wrap(~ Metric, scales = "free_x", nrow = 3) +
  labs(x = "Metric value",
       y = NULL) +
  theme_ridges(grid = TRUE, center_axis_labels = TRUE) +
  theme(legend.position    = "none",
        strip.text         = element_text(face = "bold", size = 9.5,
                                          colour = "white"),
        strip.background   = element_rect(fill = "grey25", colour = NA),
        axis.text.y        = element_text(face = "italic", size = 9),
        axis.text.x        = element_text(size = 7.5),
        axis.title.x       = element_text(face = "bold", size = 10,
                                          margin = margin(t = 6)),
        panel.spacing.x    = unit(1.0, "lines"),
        panel.spacing.y    = unit(0.8, "lines"),
        plot.margin        = margin(6, 8, 6, 6))

save_publication_plot(plot_ridge_conf,
                      file_name = "Figure3_Ridgeline_Configuration.tiff",
                      width_mm  = 174, height_mm = 200)


# ==============================================================================
# 9. EXPORT MEDIAN SUMMARY (for figure captions and Results text)     [7/9]
# ==============================================================================
message("\n[7/9] Exporting median summary tables (for caption and Results)...")

median_comp_wide <- median_comp %>%
  pivot_wider(names_from = Metric, values_from = Median)
write.csv(median_comp_wide,
          file.path(output_path, "Table3_Medians_Composition_PLAND.csv"),
          row.names = FALSE)

median_conf_wide <- df_long_conf %>%
  group_by(Metric, GENUS) %>%
  summarise(Median = median(Value, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Metric, values_from = Median)
write.csv(median_conf_wide,
          file.path(output_path, "Table4_Medians_Configuration.csv"),
          row.names = FALSE)

# ------------------------------------------------------------------------------
# TRACEABILITY OF THE SUPPLEMENTARY DESCRIPTIVE TABLES (added 7 Aug 2026)
#
# Supplementary Tables S3, S4 and S6 carry three blocks of statistics that no
# script exported: the median Percentage of Landscape of the Non-Vegetated class
# by genus, the mean, standard deviation and median of the Herbaceous and
# Non-Vegetated classes by species, and the mean and median Largest Patch Index
# of forest by species. They were correct, having been computed from the same
# file, but they could not be traced to an exported table, which is a
# reproducibility gap in a manuscript that claims code availability. They are
# exported here so that every number in the supplement has a file behind it.
# ------------------------------------------------------------------------------
write.csv(
  data_raw %>%
    group_by(GENUS) %>%
    summarise(PLAND_NonVegetated_Median = median(PLAND_Non_Vegetated, na.rm = TRUE),
              PLAND_Water_Median        = median(PLAND_Water, na.rm = TRUE),
              .groups = "drop"),
  file.path(output_path, "Table3b_Medians_MinorClasses.csv"), row.names = FALSE)

write.csv(
  data_raw %>%
    group_by(SPECIES) %>%
    summarise(
      N_SUs                     = n(),
      PLAND_Herbaceous_Mean     = mean(PLAND_Herbaceous, na.rm = TRUE),
      PLAND_Herbaceous_SD       = sd(PLAND_Herbaceous, na.rm = TRUE),
      PLAND_Herbaceous_Median   = median(PLAND_Herbaceous, na.rm = TRUE),
      PLAND_NonVegetated_Mean   = mean(PLAND_Non_Vegetated, na.rm = TRUE),
      PLAND_NonVegetated_SD     = sd(PLAND_Non_Vegetated, na.rm = TRUE),
      PLAND_NonVegetated_Median = median(PLAND_Non_Vegetated, na.rm = TRUE),
      LPI_Forest_Mean           = mean(LPI_Forest, na.rm = TRUE),
      LPI_Forest_Median         = median(LPI_Forest, na.rm = TRUE),
      PLAND_Forest_Min          = min(PLAND_Forest, na.rm = TRUE),
      PLAND_Forest_Max          = max(PLAND_Forest, na.rm = TRUE),
      .groups = "drop"),
  file.path(output_path, "Table1b_MinorClasses_And_LPI_By_Species.csv"), row.names = FALSE)

message("   Table3b and Table1b exported: every supplementary descriptive value now has a source file.")


# ==============================================================================
# 9.5 SAMPLE SIZE ACTUALLY USED BY EACH ANALYSIS OF PART I
# ==============================================================================
# The number of units entering each analysis is a reportable quantity, not an
# implementation detail. It is written to disk so that the Methods, the figure
# captions and the supplementary tables can quote it from a single source
# instead of being kept in step by hand.
# ==============================================================================
sample_sizes_I <- data.frame(
  Analysis = c("Rows loaded",
               "NMDS on Bray-Curtis of sqrt(CA)",
               "PERMANOVA and betadisper",
               "envfit of configuration metrics",
               "Units with centroid coordinates",
               "Complete over the analysis pool"),
  N        = c(nrow(data_raw),
               nrow(LULC_sqrt),
               nrow(as.matrix(LULC_dist)),
               n_envfit,
               sum(!is.na(data_raw$coord_x)),
               sum(complete.cases(data_raw[, pool_present])))
)
write.csv(sample_sizes_I,
          file.path(output_path, "TableS_Sample_Sizes_PartI.csv"),
          row.names = FALSE)
message("\n   Sample size per analysis:")
print(sample_sizes_I, row.names = FALSE)

if (any(sample_sizes_I$N != N_UNITS_EXPECTED)) {
  stop("At least one Part I analysis did not use all ", N_UNITS_EXPECTED,
       " sampling units. See the table above.")
}


# ==============================================================================
# 10. EXPORT DATA WITH COORDINATES for Part II                         [8/9]
# ==============================================================================
# Part II requires Data_Raw_WithCoords.csv and will stop if it is missing.
# The integrity check below confirms that coord_x / coord_y were populated
# by Section 4 before writing.
# ==============================================================================
message("\n[8/9] Exporting Data_Raw_WithCoords.csv...")

if (!all(c("coord_x", "coord_y") %in% names(data_raw))) {
  stop(
    "coord_x / coord_y columns are missing from data_raw.\n",
    "Check Section 4 (centroid extraction) for errors before proceeding."
  )
}

pct_with_coords <- 100 * mean(!is.na(data_raw$coord_x))
if (pct_with_coords < 100) {
  print(as.data.frame(
    data_raw[is.na(data_raw$coord_x), c("SPECIES", "UA_ID")]), row.names = FALSE)
  stop(sprintf(paste0("Only %.1f%% of rows have coordinates. Part II would drop ",
                      "the rest from the MEM screening, the variation ",
                      "partitioning and the RDA. Resolve them in Section 4."),
               pct_with_coords))
}

# Final guard before the file leaves Part I. Part II and Part III both take this
# file as their sole input, so it must already satisfy the sample-size contract.
stopifnot(
  nrow(data_raw) == N_UNITS_EXPECTED,
  sum(complete.cases(data_raw[, pool_present])) == N_UNITS_EXPECTED,
  all(!is.na(data_raw$coord_x)), all(!is.na(data_raw$coord_y))
)

write.csv(data_raw,
          file.path(data_path, "Data_Raw_WithCoords.csv"),
          row.names = FALSE)
message(sprintf(
  "   Data_Raw_WithCoords.csv exported (%d rows, %.1f%% with coordinates).",
  nrow(data_raw), pct_with_coords
))


# ==============================================================================
# 11. SESSION INFO                                                      [9/9]
# ==============================================================================
message("\n[9/9] All Part I analyses completed.")
message("\nOutputs in: ", output_path)
message("  Figures:")
message("    Figure1_NMDS_Composition_Genus.tiff")
message("    Figure2_Ridgeline_Composition.tiff      (Forest / Farming / Shrub/Herb)")
message("    Figure3_Ridgeline_Configuration.tiff")
message("    FigS2_Composition_Correlation.tiff      (Spearman, 3 focal classes)")
message("  Tables:")
message("    Table1_Composition_By_Species.csv        (PLAND_Forest, PLAND_Farming, PLAND_ShrubHerb)")
message("    Table2_Configuration_By_Species.csv      (LPI removed)")
message("    Table3_Medians_Composition_PLAND.csv")
message("    Table4_Medians_Configuration.csv")
message("    TableS2_Envfit_Configuration.csv")
message("    TableS3_Betadisper.csv")
message("    TableS4_PERMANOVA_Global.csv")
message("    TableS5_PERMANOVA_Pairwise_{Genus,Species,Diet,Locomotion}.csv")
message("    TableS_Metric_Repair_Log.csv            (cells filled from Part 0)")
message("    TableS_Sample_Sizes_PartI.csv           (n per analysis; all = 67)")
message("\nData forwarded to Part II:")
message("    ", file.path(data_path, "Data_Raw_WithCoords.csv"))
message("    (coord_x / coord_y populated in Section 4; required by Part II)")
message("\nSample size: ", N_UNITS_EXPECTED,
        " sampling units in every analysis of Part I, and in the file")
message("forwarded to Parts II and III. Taxon labels are canonical; note that")
message("Brachyteles hypoxanthus replaces the earlier spelling B_hypoxantus.")

message("\n[Session]")
message("\n[Session] Package versions and system information:")
print(sessionInfo())
writeLines(capture.output(sessionInfo()),
           file.path(output_path, "session_info_part01.txt"))

# ==============================================================================
# END OF PART I SCRIPT (v10)
# ==============================================================================
