# ==============================================================================
# 02  Constrained ordination and variation partitioning
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
#   Identify which landscape and environmental predictors structure composition,
#   and how much of the variation in envelope extent is attributable to the
#   landscape once spatial structure is accounted for. Redundancy analysis is
#   run on the Hellinger-transformed class area of the three focal classes, and
#   the term tests it reports are marginal. Variation partitioning of
#   log(extent) separates a landscape set from a spatial set built from
#   distance-based Moran eigenvector maps; individual eigenvectors are selected
#   only if the global spatial model is significant, which controls the Type I
#   error inflation that arises otherwise. Collinearity is reduced iteratively at
#   a variance inflation threshold of 5, following Zuur et al. (2010).
#
# REVISION OF 22 AUGUST 2026
#   Three changes, each documented at the point where it applies.
#   1. Section 11. The ordination response was the square root of the ABSOLUTE
#      class area, which carries the area of the sampling unit and therefore the
#      response variable of Part III: its first principal component accounted for
#      72.5% of the total inertia and correlated at r = 0.993 with the square
#      root of the unit area. It is now Hellinger-transformed, which removes the
#      size component (Legendre & Gallagher 2001).
#   2. Section 11. Term tests were sequential, by = "terms", and were reported as
#      marginal. Marginal tests are now computed and reported; the sequential
#      decomposition is exported beside them.
#   3. Section 9. The fractions of varpart() were read by position, which
#      exchanged the shared and the purely spatial fractions. They are now read
#      by the row names vegan writes, and two arithmetic identities are asserted
#      before anything is exported.
#
# INPUT   Dados/Processados/Data_Raw_WithCoords.csv
# OUTPUT  Outputs/Manuscrito/PartII/, including TableII_2b_MEM_Vectors.csv,
#         which script 03 reads for its spatial sensitivity check
#
# GUARDS
#   Verifies the input rather than repairing it. Stops if any unit is incomplete
#   over the analysis pool, if any lacks coordinates, or if edge density is zero
#   in a unit where the class is present.
#
# FIGURE TITLES
#   Titles belong in the captions, not inside the images. See script 01.
# ==============================================================================


# ==============================================================================
# 1. ENVIRONMENT SETUP
# ==============================================================================
suppressPackageStartupMessages({
  library(here)          # reproducible relative paths
  library(sf)            # spatial vector data (safeguard join only)
  library(tidyverse)     # data wrangling and ggplot2
  library(vegan)         # varpart, rda, anova.cca, RsquareAdj
  library(adespatial)    # dbmem, forward.sel
  library(car)           # vif
  library(corrplot)      # correlation plots
  library(ggplot2)       # graphics
  library(ggrepel)       # non-overlapping labels in RDA biplot
})

# Reproducibility
GLOBAL_SEED <- 123
N_PERM      <- 9999

# The full sampling design. Part II must run on this many units.
N_UNITS_EXPECTED <- 67

# Same definition as Part 0 and Part I. Kept here so that Part II can be run on
# an input produced by an older version of Part I without silently mismatching
# taxon labels. Keep the three copies in step.
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

# Paths (relative to project root)
data_path   <- here("Dados", "Processados")
output_path <- here("Outputs", "Manuscrito", "PartII")

if (!dir.exists(output_path)) dir.create(output_path, recursive = TRUE)


# ==============================================================================
# 2. DATA LOADING                                                      [1/9]
# ==============================================================================
# Requires Data_Raw_WithCoords.csv produced by Part I (v7 or later).
# Run Part I first if this file does not exist.
# ==============================================================================
message("\n[1/9] Loading Data_Raw_WithCoords.csv (produced by Part I)...")

coords_file <- file.path(data_path, "Data_Raw_WithCoords.csv")
if (!file.exists(coords_file)) {
  stop(
    "Data_Raw_WithCoords.csv not found in: ", data_path, "\n",
    "Run Part I (v7 or later) first to generate this file."
  )
}

data_raw <- read_csv(coords_file, show_col_types = FALSE)
message(sprintf("   Loaded: %d rows.", nrow(data_raw)))

# Backwards compatibility for column rename
if ("UA_Area_ha" %in% names(data_raw) &&
    !("AOH_unit_area_ha" %in% names(data_raw))) {
  data_raw <- data_raw %>% rename(AOH_unit_area_ha = UA_Area_ha)
}

data_raw <- data_raw %>%
  mutate(across(c(ORDER, GENUS, SPECIES, DIET, LOCOMOTION), as.factor),
         UA_ID = as.integer(UA_ID))

# ------------------------------------------------------------------------------
# INPUT VERIFICATION (v6). Part II no longer repairs anything.
#
# Two operations used to live here and have moved to Part I, Section 3.5, where
# the primary data are still in one place: filling the empty cells of A_guariba
# UA 6 and B_arachnoides UA 1 from the recomputation of Part 0, and recoding
# structural zeros. The recoding is now conditional on the class area being zero.
# The previous unconditional version wrote ED_Forest = 0 into two units that hold
# more than 90 per cent forest cover, and those false zeros propagated into the
# variation partitioning and into the RDA without any diagnostic showing it.
#
# What remains here is a check. If the input does not already satisfy the
# sample-size contract, Part II stops and says which cells are at fault, instead
# of quietly analyzing a subset.
# ------------------------------------------------------------------------------
data_raw$SPECIES <- canonical_taxon(data_raw$SPECIES)
data_raw$UA_ID   <- as.integer(data_raw$UA_ID)

if (nrow(data_raw) != N_UNITS_EXPECTED) {
  stop(sprintf(paste0("Data_Raw_WithCoords.csv holds %d rows, expected %d. ",
                      "Rerun Part I (v10 or later)."),
               nrow(data_raw), N_UNITS_EXPECTED))
}

ANALYSIS_POOL <- c(
  "PLAND_Forest", "PLAND_Agropecuaria", "PLAND_Herbaceous",
  "PD_Forest", "ED_Forest", "Mean_FRAC_Forest", "Mean_ENN_Forest_m",
  "PROX_MN_Forest", "Mean_AREA_Forest_ha",
  "PD_Agropecuaria", "ED_Agropecuaria", "PROX_MN_Agropecuaria",
  "Mean_Elevation_m", "Mean_Pop_Density"
)
pool_present <- intersect(ANALYSIS_POOL, names(data_raw))
missing_pool <- setdiff(ANALYSIS_POOL, names(data_raw))
if (length(missing_pool) > 0) {
  stop("Predictor(s) absent from the input: ", paste(missing_pool, collapse = ", "))
}

incomplete <- !complete.cases(data_raw[, pool_present])
if (any(incomplete)) {
  offenders <- data_raw %>%
    filter(incomplete) %>%
    select(SPECIES, UA_ID, all_of(pool_present)) %>%
    pivot_longer(-c(SPECIES, UA_ID), names_to = "Column", values_to = "Value") %>%
    filter(is.na(Value))
  print(as.data.frame(offenders), row.names = FALSE)
  stop(sprintf(paste0("%d unit(s) are incomplete over the analysis pool and ",
                      "would be dropped by complete.cases(). Run ",
                      "Script_Part_0_Recompute_Units_v2.R and then Part I ",
                      "(v10 or later)."),
               sum(incomplete)))
}

if (!all(c("coord_x", "coord_y") %in% names(data_raw)) ||
    any(is.na(data_raw$coord_x)) || any(is.na(data_raw$coord_y))) {
  stop("Some units lack centroid coordinates. Part I (v10 or later) stops on ",
       "this condition; rerun it rather than filtering here.")
}

# A false zero is harder to see than a missing value, so the two variables that
# carried them are checked explicitly against the class area.
for (cls in c("Forest", "Agropecuaria", "Herbaceous")) {
  ca_col <- paste0("CA_", cls); ed_col <- paste0("ED_", cls)
  if (all(c(ca_col, ed_col) %in% names(data_raw))) {
    suspect <- data_raw[[ca_col]] > 0 & data_raw[[ed_col]] == 0
    if (any(suspect, na.rm = TRUE)) {
      print(as.data.frame(data_raw[which(suspect), c("SPECIES", "UA_ID",
                                                     ca_col, ed_col)]),
            row.names = FALSE)
      stop(sprintf(paste0("ED_%s is zero in unit(s) where the class is present. ",
                          "This is the signature of the unconditional structural ",
                          "zero removed in v6. Rerun Part I (v10 or later)."), cls))
    }
  }
}

message(sprintf("   Input verified: %d units, complete over %d predictors, all with coordinates.",
                nrow(data_raw), length(pool_present)))
message(sprintf("   %d observations across %d species and %d genera.",
                nrow(data_raw),
                n_distinct(data_raw$SPECIES),
                n_distinct(data_raw$GENUS)))


# ==============================================================================
# 3. PEARSON CORRELATION MATRIX OF CANDIDATE PREDICTORS               [2/9]
# ==============================================================================
# Exports FigS1 before any variable selection, giving a full picture of
# pairwise collinearity among all candidates entering the VP and RDA pools.
# ==============================================================================
message("\n[2/9] Building Pearson correlation matrix of candidate predictors...")

all_candidates <- c(
  "PLAND_Forest", "PLAND_Agropecuaria", "PLAND_Herbaceous",
  "PD_Forest",    "ED_Forest",          "Mean_FRAC_Forest",
  "Mean_ENN_Forest_m", "PROX_MN_Forest", "Mean_AREA_Forest_ha",
  "PD_Agropecuaria", "ED_Agropecuaria",  "PROX_MN_Agropecuaria",
  "Mean_Elevation_m", "Mean_Pop_Density"
)

# Keep only candidates present in data_raw
all_candidates <- intersect(all_candidates, names(data_raw))

# Section 2 already established that every unit carries coordinates, so no
# filtering happens here. The object is kept under its old name so that the rest
# of the script reads unchanged.
data_sp_full <- data_raw
n_sp <- nrow(data_sp_full)
stopifnot(n_sp == N_UNITS_EXPECTED)
message(sprintf("   %d sampling units with coordinates.", n_sp))

cor_mat_candidates <- cor(data_sp_full[, all_candidates],
                          use = "complete.obs", method = "pearson")

tiff(filename    = file.path(output_path, "FigS1_Collinearity_Pearson.tiff"),
     width       = 174, height = 174, units = "mm",
     res         = 600, compression = "lzw")
corrplot(cor_mat_candidates,
         method      = "number", type = "upper",
         tl.col      = "black",  tl.srt = 45,
         number.cex  = 0.60,
         mar         = c(0, 0, 2, 0),
         title       = "Pearson correlation matrix of landscape predictors")
invisible(dev.off())
message("   FigS1_Collinearity_Pearson.tiff exported.")


# ==============================================================================
# 4. ANALYTICAL SAMPLE SIZE
# ==============================================================================
# The emergency centroid rebuild that used to sit here has been removed. It
# existed because Part I could finish with missing coordinates and only warn; it
# duplicated the centroid logic in a second place, and the two copies drifted
# apart. Part I (v10 or later) now stops on that condition, so the fallback is
# unreachable and its removal takes one source of divergence out of the pipeline.
#
# The number of units entering each analysis is a reportable quantity. It is
# printed here and written to disk in Section 11, so that the Methods and the
# supplementary tables quote one source rather than being kept in step by hand.
# ==============================================================================
data_sp <- data_raw
n_sp    <- nrow(data_sp)

stopifnot(n_sp == N_UNITS_EXPECTED,
          all(!is.na(data_sp$coord_x)), all(!is.na(data_sp$coord_y)))
message(sprintf("\n   %d sampling units retained for spatial analyses, all with coordinates.",
                n_sp))


# ==============================================================================
# 5. PREDICTOR POOLS                                                   [3/9]
# ==============================================================================
# VP pool:  composition (PLAND) + forest config + matrix config + environment.
#           PLAND is a legitimate predictor of the univariate response log(AOH).
# RDA pool: forest config + matrix config + environment only.
#           PLAND is excluded because the RDA response IS the LULC composition
#           (using PLAND as a predictor would be tautological).
#
# Exclusions applied:
#   - LPI_Forest removed (VIF > 40; redundant with PLAND_Forest)
#   - LPI_Agropecuaria removed (same reason)
#   - PLAND_Non_Vegetated removed (class excluded at KDE stage)
#   - CA and NP excluded (covary mechanically with sampling-unit size)
#   - Water class excluded (median PLAND = 0.25%)
# ==============================================================================
message("\n[3/9] Defining predictor pools...")

# 5.1 Composition (3 focal LULC classes)
pool_composition <- c(
  "PLAND_Forest",
  "PLAND_Agropecuaria",
  "PLAND_Herbaceous"
)

# 5.2 Forest configuration (LPI_Forest removed)
pool_config_forest <- c(
  "PD_Forest",           # patch density
  "ED_Forest",           # edge density
  "Mean_FRAC_Forest",    # mean fractal dimension (shape complexity)
  "Mean_ENN_Forest_m",   # mean Euclidean nearest-neighbour distance
  "PROX_MN_Forest",      # mean proximity index
  "Mean_AREA_Forest_ha"  # mean patch area
)

# 5.3 Agropecuária / farming matrix configuration (LPI_Agropecuaria removed)
pool_config_matrix <- c(
  "PD_Agropecuaria",
  "ED_Agropecuaria",
  "PROX_MN_Agropecuaria"
)

# 5.4 Environmental covariates
pool_environment <- c(
  "Mean_Elevation_m",
  "Mean_Pop_Density"
)

# 5.5 VP pool: composition + configuration + environment
pool_VP  <- c(pool_composition, pool_config_forest, pool_config_matrix, pool_environment)

# 5.6 RDA pool: configuration + environment (no PLAND)
pool_RDA <- c(pool_config_forest, pool_config_matrix, pool_environment)

message(sprintf("   VP pool: %d predictors | RDA pool: %d predictors.",
                length(pool_VP), length(pool_RDA)))

# Sanity check
missing_pred <- setdiff(unique(c(pool_VP, pool_RDA)), names(data_sp))
if (length(missing_pred) > 0) {
  stop("Predictors missing from data: ", paste(missing_pred, collapse = ", "))
}


# ==============================================================================
# 6. COLLINEARITY DIAGNOSTICS — VIF REDUCTION                         [4/9]
# ==============================================================================
# Reference: Zuur et al. (2010) A protocol for data exploration to avoid
# common statistical problems. Methods in Ecology and Evolution 1:3-14.
#
# Three strategies are available:
#
#   STRATEGY 1 (DEFAULT — applied): iterative VIF.
#     The predictor with the highest VIF is dropped and the model is refitted
#     until all remaining VIFs are <= threshold (5). LPI was already removed
#     from the pools above (manual step; Section 5), which is the first
#     substep of Strategy 1.
#
# Two alternatives were carried in the file as inert code and were removed on
# 24 August 2026: a principal component analysis of the predictor pool, and a
# forward selection applied after the variance inflation factor. Neither ever
# ran. The forward selection of Blanchet et al. (2008) is used in this script,
# but on the Moran eigenvectors of Section 8, which is where it is needed.
# ==============================================================================
message("\n[4/9] Iterative VIF reduction...")

vif_iterative <- function(df, vars, response_expr,
                          threshold = 5, max_iter = 30) {
  retained <- vars
  for (it in seq_len(max_iter)) {
    f <- as.formula(
      paste0(response_expr, " ~ ", paste(retained, collapse = " + "))
    )
    m <- tryCatch(lm(f, data = df), error = function(e) NULL)
    if (is.null(m)) {
      stop("VIF model failed to fit. Check for perfectly collinear inputs.")
    }
    v <- tryCatch(car::vif(m), error = function(e) NULL)
    if (is.null(v) || any(!is.finite(v))) {
      cor_mat  <- cor(df[, retained], use = "complete.obs")
      drop_var <- retained[which.max(rowSums(abs(cor_mat)) - 1)]
      message(sprintf("   Iter %02d: perfect collinearity; dropping %s", it, drop_var))
    } else {
      if (max(v) <= threshold) break
      drop_var <- names(v)[which.max(v)]
      message(sprintf("   Iter %02d: dropping %s (VIF = %.2f)", it, drop_var, max(v)))
    }
    retained <- setdiff(retained, drop_var)
    if (length(retained) < 2) break
  }
  # ADDED 24 Aug 2026. v holds the factors of the model of the previous iteration.
  # When the loop ends because fewer than two predictors remain, or because
  # max_iter was reached, that model is not the retained set, and the stale vector
  # used to be written straight into TableII_1_VIF_*.csv. Recomputing it here
  # guarantees that the exported factors describe the set that was kept.
  v_final <- tryCatch(
    car::vif(lm(as.formula(paste(response_expr, "~",
                                 paste(retained, collapse = " + "))), data = df)),
    error = function(e) NULL)
  list(retained = retained, vif_final = if (is.null(v_final)) v else v_final)
}

# 6.1 Reduce VP pool
message("\n   VP pool (response: log(AOH_unit_area_ha))")
vif_VP        <- vif_iterative(data_sp, pool_VP,
                               response_expr = "log(AOH_unit_area_ha)")
predictors_VP <- vif_VP$retained

# 6.2 Reduce RDA pool (proxy: PCoA axis 1 of Bray-Curtis composition)
message("\n   RDA pool (response proxy: PCoA1 of compositional dissimilarity)")
comp_mat_vif <- data_sp %>%
  transmute(Forest     = sqrt(replace_na(CA_Forest,       0)),
            Farming    = sqrt(replace_na(CA_Agropecuaria, 0)),
            ShrubHerb  = sqrt(replace_na(CA_Herbaceous,   0))) %>%
  as.matrix()
bray_dist  <- vegan::vegdist(comp_mat_vif, method = "bray")
pcoa1      <- cmdscale(bray_dist, k = 1)[, 1]
data_sp$comp_pcoa1 <- pcoa1

vif_RDA        <- vif_iterative(data_sp, pool_RDA,
                                response_expr = "comp_pcoa1")
predictors_RDA <- vif_RDA$retained

# Export VIF tables
write.csv(
  data.frame(Variable = names(vif_VP$vif_final),
             VIF      = round(as.numeric(vif_VP$vif_final), 3),
             Pool     = "Variation Partitioning"),
  file.path(output_path, "TableII_1_VIF_VP.csv"), row.names = FALSE
)
write.csv(
  data.frame(Variable = names(vif_RDA$vif_final),
             VIF      = round(as.numeric(vif_RDA$vif_final), 3),
             Pool     = "RDA constrained ordination"),
  file.path(output_path, "TableII_1_VIF_RDA.csv"), row.names = FALSE
)

message(sprintf("\n   Predictors retained — VP (n = %d): %s",
                length(predictors_VP), paste(predictors_VP, collapse = ", ")))
message(sprintf("   Predictors retained — RDA (n = %d): %s",
                length(predictors_RDA), paste(predictors_RDA, collapse = ", ")))


# ------------------------------------------------------------------------------
# METHODOLOGICAL NOTE: PD_Forest x Mean_FRAC_Forest
# ------------------------------------------------------------------------------
# The correlation between patch density and mean fractal dimension of the forest
# class exceeds the conventional 0.7 screening threshold (Dormann et al. 2013).
# Both predictors are nevertheless retained, because their variance inflation
# factors remain below 5 in the multivariable model. Following Zuur et al. (2010),
# VIF is the formal multivariable diagnostic; a pairwise correlation only flags a
# pair for inspection.
#
# The biological reading adopted in the Discussion treats the two jointly as one
# morphological integrity signature of the forest remnant. Landscapes with more,
# smaller and more irregularly shaped patches share a single underlying process,
# fragmentation by anthropological conversion, and the two metrics quantify
# different facets of it. They are not interpreted as independent drivers.
#
# The coefficients are computed here rather than written into the comment. The
# values previously recorded in this file, r = 0.659 and rho = 0.723, were
# obtained on 65 units, before the two units of Section 2 re-entered the
# analysis. Quote the printed values, and the printed n, in Section 2.5.5 of the
# manuscript and in the legend of Supplementary Table S21.
# ------------------------------------------------------------------------------
cor_pd_frac <- data_sp %>%
  select(PD_Forest, Mean_FRAC_Forest) %>%
  filter(complete.cases(.))

pearson_pd_frac  <- cor(cor_pd_frac$PD_Forest, cor_pd_frac$Mean_FRAC_Forest,
                        method = "pearson")
spearman_pd_frac <- cor(cor_pd_frac$PD_Forest, cor_pd_frac$Mean_FRAC_Forest,
                        method = "spearman")

message(sprintf("\n   PD_Forest x Mean_FRAC_Forest on n = %d: Pearson r = %.3f | Spearman rho = %.3f",
                nrow(cor_pd_frac), pearson_pd_frac, spearman_pd_frac))

write.csv(
  data.frame(Pair    = "PD_Forest x Mean_FRAC_Forest",
             N       = nrow(cor_pd_frac),
             Pearson = round(pearson_pd_frac, 4),
             Spearman = round(spearman_pd_frac, 4)),
  file.path(output_path, "TableII_1c_Correlation_PD_FRAC.csv"), row.names = FALSE
)


# Two optional reduction strategies, a principal component analysis and a
# forward selection, were carried here as inert code until 24 August 2026. The
# first was commented out and the second sat inside if (FALSE) and referred to
# data_z, which Section 7 creates below it, so enabling it would have failed.
# Neither ever ran, neither produced any value the manuscript reports, and the
# file TableII_1b_ForwardSel_VP.csv the second one names was never written. The
# reduction the analysis actually uses is the iterative variance inflation
# factor of Section 6, and the forward selection of Blanchet et al. (2008) is
# applied where it belongs, to the Moran eigenvectors in Section 8.


# ==============================================================================
# 7. STANDARDIZATION
# ==============================================================================
data_z          <- data_sp
vars_to_scale   <- unique(c(predictors_VP, predictors_RDA))
data_z[, vars_to_scale] <- scale(data_sp[, vars_to_scale])


# ==============================================================================
# 8. SPATIAL EIGENVECTORS (MEMs) — Moran Eigenvector Maps             [5/9]
# ==============================================================================
# WHAT THIS BLOCK TESTS, AND WHAT IT DOES NOT (added 29 Jul 2026)
#
# Scope. dbmem() builds distance-based Moran Eigenvector Maps from the centroid
# coordinates of the sampling units. The global model below regresses log(AOH)
# on the FULL set of eigenvectors and tests it by permutation. This screens for
# spatial autocorrelation AMONG sampling units, i.e. the tendency of units that
# are geographically close to hold similar values of the response.
#
# It does NOT test grouping dependence among units that share the same landscape
# unit (UA_ID). That is a different form of non-independence and is screened
# separately in Part III, section 6.5, by comparing the OLS model with a linear
# mixed model carrying a random intercept for UA_ID (ICC + AICc).
#
# The two diagnostics must be reported as answering two distinct questions.
# Presenting them as independent confirmations of the same fact is inaccurate.
#
# Positive autocorrelation only. MEM.autocor = "positive" generates eigenvectors
# modelling positive autocorrelation exclusively. This is the form that inflates
# Type I error rates in this design; negative autocorrelation is not screened,
# and the manuscript should state this explicitly.
#
# Two-step selection. Forward selection of individual eigenvectors runs only if
# the global test is significant (Blanchet et al. 2008). This controls the Type I
# error inflation that arises when eigenvectors are selected without a prior
# global test. When the global model is not significant, no eigenvector is
# retained and the spatial fraction of the partition is null by construction.
#
# FOR THE MANUSCRIPT: mem_global_p and mem_global_R2a (computed just below) are
# the two numbers the reviewer asked for and that Results par. 31 does not yet
# report. Print them and paste them into the text.
# ==============================================================================
message("\n[5/9] Building MEMs with adespatial::dbmem and forward selection...")

y_log      <- log(data_z$AOH_unit_area_ha)
coords_xy  <- as.matrix(data_z[, c("coord_x", "coord_y")])

set.seed(GLOBAL_SEED)
mem_full <- dbmem(coords_xy, MEM.autocor = "positive", silent = TRUE)
message(sprintf("   %d positive MEMs generated by dbmem.", ncol(mem_full)))

mem_global     <- rda(y_log ~ ., data = as.data.frame(mem_full))
set.seed(GLOBAL_SEED)
mem_global_anova <- anova(mem_global, permutations = N_PERM)
mem_global_p     <- mem_global_anova$`Pr(>F)`[1]
mem_global_R2a   <- RsquareAdj(mem_global)$adj.r.squared

if (is.na(mem_global_p) || mem_global_p > 0.05) {
  message(sprintf(
    "   Global spatial model not significant (p = %.4f). No MEM retained.",
    mem_global_p
  ))
  mems_selected <- NULL
  mem_mat       <- matrix(numeric(0), nrow = n_sp, ncol = 0)
} else {
  message(sprintf(
    "   Global spatial model significant (p = %.4f, R2a = %.3f).",
    mem_global_p, mem_global_R2a
  ))
  set.seed(GLOBAL_SEED)
  mem_sel <- forward.sel(y_log, as.matrix(mem_full),
                         alpha    = 0.05,
                         R2thresh = max(mem_global_R2a, 0.001),
                         nperm    = N_PERM)
  mems_selected <- mem_sel$variables
  mem_mat       <- as.matrix(mem_full[, mems_selected, drop = FALSE])
  message(sprintf("   MEMs retained: %d (%s)",
                  length(mems_selected), paste(mems_selected, collapse = ", ")))
  write.csv(mem_sel,
            file.path(output_path, "TableII_2_MEMs_Selected.csv"),
            row.names = FALSE)

  # ----------------------------------------------------------------------------
  # The eigenvectors themselves are exported, not only their names. Script 03
  # needs them to refit the confirmatory model with the spatial structure as a
  # covariate, which is the sensitivity check that Diniz-Filho, Rangel and Bini
  # (2008) require whenever model selection is carried out on data with detected
  # spatial autocorrelation: an information criterion computed on a non-spatial
  # regression produces unstable and overfitted minimum adequate models.
  # The key (SPECIES, UA_ID) is carried so that the join in script 03 cannot
  # rely on row order.
  # ----------------------------------------------------------------------------
  write.csv(
    data.frame(SPECIES = as.character(data_sp$SPECIES),
               UA_ID   = data_sp$UA_ID,
               mem_mat, check.names = FALSE),
    file.path(output_path, "TableII_2b_MEM_Vectors.csv"), row.names = FALSE)
  message("   TableII_2b_MEM_Vectors.csv exported for the spatial sensitivity check of script 03.")
}


# ==============================================================================
# 9. VARIATION PARTITIONING (adjusted R-squared)                      [6/9]
# ==============================================================================
# The fractions of varpart() are read by ROW NAME, not by position. Reading them
# by position was the defect corrected on 22 August 2026: in vegan 2.6.4 row 2 of
# $part$indfract is the purely spatial fraction and row 3 is the shared one, so
# positional reading exchanged the two. Section 9.2 asserts the two arithmetic
# identities that any two-set partition must satisfy before anything is written.
# Permutation tests use the direct rda(Y, X, Z) interface.
# ==============================================================================
message("\n[6/9] Variation Partitioning: landscape (X1) vs space (X2)...")

# Filter to complete cases for VP predictors before any rda()/varpart() call.
# Configuration metrics such as PROX_MN and ENN_MN can be NA when a class
# has fewer than two patches in a given sampling unit; rda() rejects NAs.
# The filter is retained as a guard, not as a working step. Before v6 it silently
# removed the two units whose PROX_MN_Agropecuaria cell was empty, which is where
# the n = 63 of the supplementary material came from. Those cells are now filled
# in Part I, so the filter must remove nothing; if it does, the run stops.
X1_all <- as.matrix(data_z[, predictors_VP])
vp_ok  <- complete.cases(X1_all)
n_excl <- sum(!vp_ok)
if (n_excl > 0) {
  print(as.data.frame(data_z[!vp_ok, c("SPECIES", "UA_ID")]), row.names = FALSE)
  stop(sprintf(paste0("%d sampling unit(s) carry an NA among the retained VP ",
                      "predictors and would be dropped. Run Part 0 and Part I ",
                      "(v10 or later) before Part II."), n_excl))
}
stopifnot(sum(vp_ok) == N_UNITS_EXPECTED)
message(sprintf("   Variation partitioning fitted on %d sampling units.", sum(vp_ok)))
y_log_vp <- y_log[vp_ok]
X1       <- X1_all[vp_ok, , drop = FALSE]
X2       <- if (ncol(mem_mat) > 0) mem_mat[vp_ok, , drop = FALSE] else mem_mat

# varpart() requires >= 2 explanatory tables. When no MEMs were retained
# (X2 is empty), the full partitioning is not possible; the landscape R-squared
# is obtained directly from rda() + RsquareAdj(), and vp_res is set to NULL
# so that Section 10 skips the Venn diagram.

if (ncol(X2) == 0) {
  
  message("   No spatial component retained — computing landscape R-squared only.")
  message("   (varpart() requires >= 2 tables; single-set handled via rda().)")
  
  set.seed(GLOBAL_SEED)
  rda_full  <- rda(y_log_vp, X1)
  test_full <- anova(rda_full, permutations = N_PERM)
  r2a_x1    <- RsquareAdj(rda_full)$adj.r.squared
  
  vp_table <- data.frame(
    Fraction       = c("[a] Landscape (X1)", "[d] Residuals (unexplained)"),
    Adj_R2         = round(c(r2a_x1, 1 - r2a_x1), 4),
    Pct_of_total   = round(100 * c(r2a_x1, 1 - r2a_x1), 2),
    Interpretation = c("Total landscape effect (no spatial component retained)",
                       "Variance not explained by landscape")
  )
  
  vp_tests <- data.frame(
    Test  = "Full landscape model (X1)",
    F_obs = round(test_full$F[1], 4),
    p_val = test_full$`Pr(>F)`[1]
  )
  
  vp_res <- NULL   # signals Section 10 to skip the Venn diagram
  
} else {
  
  # Two-set partitioning: varpart() is safe to call
  vp_res <- varpart(y_log_vp, X1, X2)
  ind    <- vp_res$part$indfract

  # --------------------------------------------------------------------------
  # CORRECTED 22 Aug 2026. Until v6 the three testable fractions were read by
  # POSITION, ind[1], ind[2], ind[3], and labelled [a] pure landscape, [b]
  # shared, [c] pure space in that order. In vegan 2.6.4 the rows of
  # $part$indfract are named "[a] = X1|X2", "[b] = X2|X1", "[c]",
  # "[d] = Residuals": row two is PURE SPACE and row three is the SHARED
  # fraction. Positional reading therefore exchanged the shared and the purely
  # spatial fractions in every exported table, and the manuscript reported
  # 1.2% shared with 9.9% purely spatial when the values are the other way round.
  #
  # Two changes prevent it from recurring. The fractions are now selected by the
  # row name that vegan itself writes, and the two arithmetic identities that any
  # correct two-set partition must satisfy are asserted before anything is
  # written to disk:
  #
  #     [a] + [b_shared] = adjusted R2 of X1 alone
  #     [b_shared] + [c] = adjusted R2 of X2 alone
  #
  # With the old positional reading the first of these failed by 0.087, which is
  # what made the error visible.
  # --------------------------------------------------------------------------
  rn         <- rownames(ind)
  row_pure1  <- grep("X1\\|X2", rn)
  row_pure2  <- grep("X2\\|X1", rn)
  row_shared <- setdiff(grep("^\\[[abc]\\]", rn), c(row_pure1, row_pure2))
  if (length(row_pure1) != 1 || length(row_pure2) != 1 || length(row_shared) != 1)
    stop("varpart() returned fraction labels this script does not recognise: ",
         paste(rn, collapse = " | "))

  frac_a      <- ind[row_pure1,  "Adj.R.squared"]   # landscape, space removed
  frac_c      <- ind[row_pure2,  "Adj.R.squared"]   # space, landscape removed
  frac_shared <- ind[row_shared, "Adj.R.squared"]   # jointly explained
  frac_d      <- 1 - (frac_a + frac_shared + frac_c)

  r2a_X1 <- RsquareAdj(rda(y_log_vp, X1))$adj.r.squared
  r2a_X2 <- RsquareAdj(rda(y_log_vp, X2))$adj.r.squared
  tol    <- 1e-6
  if (abs((frac_a + frac_shared) - r2a_X1) > tol ||
      abs((frac_shared + frac_c) - r2a_X2) > tol) {
    stop(sprintf(paste0("Variation partitioning failed its identity check.\n",
                        "  [a] + shared = %.6f, adjusted R2 of X1 alone = %.6f\n",
                        "  shared + [c] = %.6f, adjusted R2 of X2 alone = %.6f\n",
                        "Inspect the row names of varpart()$part$indfract."),
                 frac_a + frac_shared, r2a_X1, frac_shared + frac_c, r2a_X2))
  }
  message(sprintf(paste0("   Partition identity check passed: [a] + shared = ",
                         "%.4f = adjR2(X1); shared + [c] = %.4f = adjR2(X2)."),
                  r2a_X1, r2a_X2))

  vp_table <- data.frame(
    Fraction       = c("[a] Pure landscape (X1 | X2)",
                       "[b] Shared landscape and space",
                       "[c] Pure space (X2 | X1)",
                       "[d] Residuals (unexplained)"),
    Adj_R2         = round(c(frac_a, frac_shared, frac_c, frac_d), 4),
    Pct_of_total   = round(100 * c(frac_a, frac_shared, frac_c, frac_d), 2),
    Interpretation = c("Landscape effect independent of geography",
                       "Landscape variation that is spatially structured",
                       "Spatial structure without landscape correspondence",
                       "Variance not explained by either set")
  )

  # Permutation tests — direct rda(Y, X, Z) interface
  set.seed(GLOBAL_SEED)
  rda_full            <- rda(y_log_vp, cbind(X1, X2))
  test_full           <- anova(rda_full, permutations = N_PERM)
  
  set.seed(GLOBAL_SEED)
  rda_pure_landscape  <- rda(y_log_vp, X1, X2)   # [a]: X1 | X2
  test_pure_landscape <- anova(rda_pure_landscape, permutations = N_PERM)
  
  set.seed(GLOBAL_SEED)
  rda_pure_space      <- rda(y_log_vp, X2, X1)   # [c]: X2 | X1
  test_pure_space     <- anova(rda_pure_space, permutations = N_PERM)
  
  vp_tests <- data.frame(
    Test  = c("Full model (X1 + X2)",
              "Pure landscape (X1 | X2)",
              "Pure space (X2 | X1)"),
    F_obs = round(c(test_full$F[1],
                    test_pure_landscape$F[1],
                    test_pure_space$F[1]), 4),
    p_val = c(test_full$`Pr(>F)`[1],
              test_pure_landscape$`Pr(>F)`[1],
              test_pure_space$`Pr(>F)`[1])
  )
}

write.csv(vp_table,
          file.path(output_path, "TableII_3_VP_Fractions.csv"),
          row.names = FALSE)
message("   VP fractions:")
print(vp_table[, c("Fraction", "Adj_R2", "Pct_of_total")], row.names = FALSE)

write.csv(vp_tests,
          file.path(output_path, "TableII_3_VP_Significance.csv"),
          row.names = FALSE)


message("   Permutation test results:")
print(vp_tests, row.names = FALSE)


# ------------------------------------------------------------------------------
# 9.1 SENSITIVITY: the partition without the one predictor tied to the response
# ------------------------------------------------------------------------------
# ADDED 22 Aug 2026. The landscape set retained after the VIF screening still
# contains PD_Forest, and Part V shows that log(OHE) = log(NP) + log(100) -
# log(PD) is an identity on these units. A fraction of variation in log(OHE)
# attributed to a set that contains PD is therefore partly arithmetic. The same
# partition is run without it, so that the manuscript can report how much of the
# landscape fraction depends on that single predictor rather than assert that it
# does not matter. Nothing here replaces the main partition; it qualifies it.
# ------------------------------------------------------------------------------
if (!is.null(vp_res) && "PD_Forest" %in% colnames(X1)) {

  X1_free <- X1[, setdiff(colnames(X1), "PD_Forest"), drop = FALSE]
  vp_free <- varpart(y_log_vp, X1_free, X2)
  indf    <- vp_free$part$indfract
  rnf     <- rownames(indf)

  af <- indf[grep("X1\\|X2", rnf), "Adj.R.squared"]
  cf <- indf[grep("X2\\|X1", rnf), "Adj.R.squared"]
  bf <- indf[setdiff(grep("^\\[[abc]\\]", rnf),
                     c(grep("X1\\|X2", rnf), grep("X2\\|X1", rnf))),
             "Adj.R.squared"]

  set.seed(GLOBAL_SEED)
  tf_full  <- anova(rda(y_log_vp, cbind(X1_free, X2)), permutations = N_PERM)
  set.seed(GLOBAL_SEED)
  tf_land  <- anova(rda(y_log_vp, X1_free, X2),        permutations = N_PERM)
  set.seed(GLOBAL_SEED)
  tf_space <- anova(rda(y_log_vp, X2, X1_free),        permutations = N_PERM)

  vp_free_tab <- data.frame(
    Fraction = c("[a] Pure landscape (X1 | X2)",
                 "[b] Shared landscape and space",
                 "[c] Pure space (X2 | X1)",
                 "[d] Residuals (unexplained)",
                 "Full model (X1 + X2)"),
    Adj_R2   = round(c(af, bf, cf, 1 - (af + bf + cf), af + bf + cf), 4),
    F_obs    = c(round(tf_land$F[1], 4), NA, round(tf_space$F[1], 4), NA,
                 round(tf_full$F[1], 4)),
    p_val    = c(tf_land$`Pr(>F)`[1], NA, tf_space$`Pr(>F)`[1], NA,
                 tf_full$`Pr(>F)`[1]),
    Note     = "landscape set without PD_Forest, the only retained metric normalised by the response"
  )
  write.csv(vp_free_tab,
            file.path(output_path, "TableII_3b_VP_Fractions_Without_PD.csv"),
            row.names = FALSE)
  message("   Partition without PD_Forest:")
  print(vp_free_tab[, c("Fraction", "Adj_R2", "F_obs", "p_val")], row.names = FALSE)
}


# ==============================================================================
# 10. VENN DIAGRAM (Figure VP)                                         [7/9]
# ==============================================================================
message("\n[7/9] Building the variation partitioning figure...")

if (is.null(vp_res)) {
  # Single-set: plot(varpart) is unavailable. Export a simple text summary.
  message("   No spatial component — Venn diagram not applicable.")
  message("   Exporting landscape-only summary figure.")
  
  tiff(filename    = file.path(output_path, "Figure_VP_Summary.tiff"),
       width       = 120, height = 80, units = "mm",
       res         = 600, compression = "lzw")
  par(mar = c(1, 1, 2.5, 1), bg = "white")
  plot.new()
  rect(0.15, 0.25, 0.85, 0.80, col = "#A6CEE3", border = "grey40", lwd = 1.2)
  text(0.50, 0.52,
       sprintf("Landscape (X1)\nAdj. R\u00B2 = %.3f",
               RsquareAdj(rda_full)$adj.r.squared),
       cex = 1.0, font = 2, adj = 0.5)
  text(0.50, 0.12,
       sprintf("Permutation p = %.4f  |  No spatial eigenvectors retained (MEM p > 0.05)",
               test_full$`Pr(>F)`[1]),
       cex = 0.72, col = "grey35", adj = 0.5)
  title("log(AOH_ha) - landscape-only model", cex.main = 0.88)
  invisible(dev.off())
  message("   Figure_VP_Summary.tiff exported.")
  
} else {
  # Two-set: standard Venn diagram via plot.varpart
  tiff(filename    = file.path(output_path, "Figure_VP_Venn.tiff"),
       width       = 174, height = 130, units = "mm",
       res         = 600, compression = "lzw")
  par(mar = c(2, 2, 3, 2))
  plot(vp_res,
       digits  = 2,
       bg      = c("#A6CEE3", "#FB9A99"),
       Xnames  = c("Landscape", "Space"),
       id.size = 1.0,
       cex     = 1.0,
       main    = "Variation partitioning of log(AOH_ha) by landscape and space")
  mtext(sprintf("Adj. R-squared (full model) = %.3f | permutation p = %.4f",
                RsquareAdj(rda_full)$adj.r.squared,
                test_full$`Pr(>F)`[1]),
        side = 1, line = 0.4, cex = 0.85)
  invisible(dev.off())
  message("   Figure_VP_Venn.tiff exported.")
}


# ==============================================================================
# 11. CONSTRAINED ORDINATION: RDA OF LULC COMPOSITION                 [8/9]
# ==============================================================================
# Response: Hellinger-transformed Class Area of the three focal LULC classes
#           (Forest, Farming, Shrub/Herb). Non-Vegetated excluded.
# Predictors: configuration and environmental variables retained after VIF
#             (PLAND excluded to avoid circularity with the response).
#
# CORRECTED 22 Aug 2026. Until v6, the response was the square root of the
# ABSOLUTE class area in hectares. Absolute class area is the product of the
# proportion of the class and the area of the sampling unit, and the area of the
# sampling unit is the response variable of Part III. The consequence is
# measurable on these data: the first principal component of sqrt(CA) accounts
# for 72.5% of the total inertia. It correlates with the square root of the unit
# area at r = 0.993. The ordination was therefore ordering the units by the size
# of the window and not by their composition, which is the same dependence that
# Part V exists to measure, transferred to the multivariate analysis.
#
# The Hellinger transformation, that is, the square root of the RELATIVE class
# area, removes the size component and is the standard preparation for a linear
# ordination of composition data (Legendre & Gallagher 2001). Under it the
# first principal component of the response correlates with the square root of
# the unit area at r = 0.054.
#
# The change also restores agreement with Part I. The unconstrained ordination
# of Part I runs on Bray-Curtis dissimilarities, which are computed on relative
# composition, and there the mean proximity index is the fourth strongest
# correlate of the composition gradient, not the first. The sqrt(CA) RDA was
# the only analysis in the chain that placed it first.
#
# Both responses are shown below. The Hellinger version is the one that
# manuscript reports; the sqrt(CA) version is kept and exported so that the
# difference between the two is comparable rather than asserted.
# ==============================================================================
message("\n[8/9] Constrained ordination: RDA of LULC composition...")

LULC_CA <- data_z %>%
  transmute(Forest       = replace_na(CA_Forest,       0),
            Farming      = replace_na(CA_Agropecuaria, 0),
            `Shrub/Herb` = replace_na(CA_Herbaceous,   0)) %>%
  as.matrix()

# Hellinger: sqrt(x_ij / row total). decostand() is used rather than a hand
# rolled formula so that the transformation is the documented one.
LULC_hell <- vegan::decostand(LULC_CA, method = "hellinger")

# Superseded response, retained for the comparison exported below.
LULC_sqrt <- sqrt(LULC_CA)

# CORRECTED 6 Aug 2026. Section 9 filtered its own copy to complete cases but
# this block did not, so rda() received NAs and stopped with
# "missing values in object". The same filter is applied here, and the number of
# units actually entering the ordination is reported rather than assumed.
X_RDA_all <- as.matrix(data_z[, predictors_RDA])
rda_ok    <- complete.cases(X_RDA_all)
if (sum(!rda_ok) > 0) {
  perdidas <- data_z %>% filter(!rda_ok) %>%
    transmute(rotulo = paste(SPECIES, UA_ID)) %>% pull(rotulo)
  stop(sprintf(paste0("%d unit(s) would be excluded from the RDA for missing ",
                      "predictors: %s. Run Part 0 and Part I (v10 or later)."),
               sum(!rda_ok), paste(perdidas, collapse = "; ")))
}
stopifnot(sum(rda_ok) == N_UNITS_EXPECTED)
message(sprintf("   RDA fitted on %d sampling units.", sum(rda_ok)))
n_rda     <- sum(rda_ok)
X_RDA     <- X_RDA_all[rda_ok, , drop = FALSE]
LULC_hell <- LULC_hell[rda_ok, , drop = FALSE]
LULC_sqrt <- LULC_sqrt[rda_ok, , drop = FALSE]

# Every object used downstream of the ordination must carry the SAME rows as the
# ordination itself. Keeping a filtered copy, rather than filtering each object
# where it is used, is what prevents the mismatch from reappearing: any block
# that reaches for data_z after this point is a bug.
data_rda <- data_z[rda_ok, , drop = FALSE]
stopifnot(nrow(data_rda) == n_rda, nrow(LULC_hell) == n_rda,
          nrow(LULC_sqrt) == n_rda, nrow(X_RDA) == n_rda)

# How much of each candidate response is the size of the window. This is the
# diagnostic that motivated the change of transformation; it is computed here so
# that the justification is a number in the output and not a claim in a comment.
pc_share <- function(Y) {
  p <- prcomp(Y, center = TRUE, scale. = FALSE)
  c(pc1_pct = 100 * p$sdev[1]^2 / sum(p$sdev^2),
    r_with_sqrt_area = abs(cor(p$x[, 1], sqrt(data_rda$AOH_unit_area_ha))))
}
size_diag <- rbind(`sqrt(class area), superseded` = pc_share(LULC_sqrt),
                   `Hellinger, reported`          = pc_share(LULC_hell))
write.csv(data.frame(Response = rownames(size_diag), round(size_diag, 4)),
          file.path(output_path, "TableII_5_Response_Size_Dependence.csv"),
          row.names = FALSE)
message(sprintf(paste0("   Size dependence of the ordination response: ",
                       "sqrt(CA) PC1 = %.1f%% of inertia, r = %.3f with ",
                       "sqrt(unit area); Hellinger PC1 = %.1f%%, r = %.3f."),
                size_diag[1, 1], size_diag[1, 2],
                size_diag[2, 1], size_diag[2, 2]))

rda_landscape    <- rda(LULC_hell ~ ., data = as.data.frame(X_RDA))

set.seed(GLOBAL_SEED)
rda_global_anova <- anova(rda_landscape, permutations = N_PERM)
# CORRECTED 22 Aug 2026. by = "terms" is a SEQUENTIAL (type I) test: each term is
# tested against what the terms before it in the formula have already explained,
# so the result depends on the order in which the predictors were written. The
# manuscript presented these tests as marginal and ordered the table by F, which
# invites the reader to compare them as if each were adjusted for all the others.
# On these data the two schemes disagree materially: under the sequential test
# forest patch density took F = 21.43 and ranked third, and under the marginal
# test on the same response it takes F = 2.96 and is not significant.
# by = "margin" is now what the manuscript reports. The sequential decomposition
# is still exported, because it is the one whose variances sum to the constrained
# inertia and it is useful for describing how much unique variation is left for
# the last terms to explain.
set.seed(GLOBAL_SEED)
rda_margin_anova <- anova(rda_landscape, by = "margin", permutations = N_PERM)
set.seed(GLOBAL_SEED)
rda_terms_anova  <- anova(rda_landscape, by = "terms", permutations = N_PERM)
set.seed(GLOBAL_SEED)
rda_axes_anova   <- anova(rda_landscape, by = "axis",  permutations = N_PERM)

rda_R2  <- RsquareAdj(rda_landscape)$r.squared
rda_R2a <- RsquareAdj(rda_landscape)$adj.r.squared

message(sprintf("   RDA: R2 = %.3f, R2-adj = %.3f, p = %.4f",
                rda_R2, rda_R2a, rda_global_anova$`Pr(>F)`[1]))

# Reported table: marginal tests, ordered by F.
rda_margin_df <- as.data.frame(rda_margin_anova) %>%
  rownames_to_column("Predictor") %>%
  arrange(desc(.data$F))
write.csv(rda_margin_df,
          file.path(output_path, "TableII_4_RDA_Anova_Marginal.csv"),
          row.names = FALSE)

rda_terms_df <- as.data.frame(rda_terms_anova) %>%
  rownames_to_column("Predictor")
write.csv(rda_terms_df,
          file.path(output_path, "TableII_4_RDA_Anova_Terms_Sequential.csv"),
          row.names = FALSE)

# Same two tests on the superseded response, exported so that the effect of the
# transformation on the ranking of the predictors can be inspected directly.
rda_superseded <- rda(LULC_sqrt ~ ., data = as.data.frame(X_RDA))
set.seed(GLOBAL_SEED)
sup_margin <- as.data.frame(anova(rda_superseded, by = "margin",
                                  permutations = N_PERM)) %>%
  rownames_to_column("Predictor")
set.seed(GLOBAL_SEED)
sup_terms  <- as.data.frame(anova(rda_superseded, by = "terms",
                                  permutations = N_PERM)) %>%
  rownames_to_column("Predictor")
write.csv(
  bind_rows(mutate(sup_margin, Test = "marginal", Response = "sqrt(class area)"),
            mutate(sup_terms,  Test = "sequential", Response = "sqrt(class area)"),
            mutate(rda_margin_df, Test = "marginal", Response = "Hellinger"),
            mutate(rda_terms_df,  Test = "sequential", Response = "Hellinger")),
  file.path(output_path, "TableII_4b_RDA_Transformation_Comparison.csv"),
  row.names = FALSE)

rda_axes_df <- as.data.frame(rda_axes_anova) %>%
  rownames_to_column("Axis")
write.csv(rda_axes_df,
          file.path(output_path, "TableII_4_RDA_Anova_Axes.csv"),
          row.names = FALSE)

# ADDED 29 Jul 2026. rda_R2, rda_R2a and the global permutation p were printed
# to the console only and never written to disk, so the
# TableII_4_RDA_GlobalSummary.csv found in Outputs was a leftover from an
# earlier run. Results paragraph 30 of the manuscript needs these three values,
# and mem_global_p / mem_global_R2a are the two numbers requested in the review.
rda_global_summary <- data.frame(
  Metric = c("R2 (constrained proportion)", "Adjusted R2",
             "Constrained inertia", "Residual inertia",
             "Global permutation p",
             "MEM global model p", "MEM global model adjusted R2"),
  Value  = c(round(rda_R2, 4), round(rda_R2a, 4),
             round(sum(rda_landscape$CCA$eig), 4),
             round(sum(rda_landscape$CA$eig), 4),
             rda_global_anova$`Pr(>F)`[1],
             mem_global_p,
             round(mem_global_R2a, 4))
)
write.csv(rda_global_summary,
          file.path(output_path, "TableII_4_RDA_GlobalSummary.csv"),
          row.names = FALSE)
message(sprintf("   RDA global exported: R2 = %.4f | adj R2 = %.4f | p = %.4f",
                rda_R2, rda_R2a, rda_global_anova$`Pr(>F)`[1]))

# Sample size actually used by each analysis, for the manuscript.
write.csv(
  data.frame(
    Analysis = c("Rows loaded", "With centroid coordinates",
                 "MEM screening", "Variation partitioning", "RDA",
                 "Correlation PD_Forest x Mean_FRAC_Forest"),
    N        = c(nrow(data_raw), n_sp, length(y_log), sum(vp_ok), n_rda,
                 nrow(cor_pd_frac))
  ),
  file.path(output_path, "TableII_0_Sample_Sizes.csv"), row.names = FALSE
)

# One number for the whole of Part II. If any analysis ran on a different sample,
# the discrepancy stops the script here rather than reaching a manuscript table.
sizes_II <- c(nrow(data_raw), n_sp, length(y_log), sum(vp_ok), n_rda,
              nrow(cor_pd_frac))
if (any(sizes_II != N_UNITS_EXPECTED)) {
  stop("Part II analyses did not all use ", N_UNITS_EXPECTED,
       " units: ", paste(sizes_II, collapse = ", "))
}
message(sprintf("   All Part II analyses ran on %d sampling units.", N_UNITS_EXPECTED))
message(sprintf("   MEM global: p = %.4f | adj R2 = %.4f",
                mem_global_p, mem_global_R2a))

# Biplot (scaling 2: relationships among predictors and response species)
rda_site_scores    <- as.data.frame(scores(rda_landscape, display = "sites",   scaling = 2))
rda_species_scores <- as.data.frame(scores(rda_landscape, display = "species", scaling = 2))
rda_biplot_arrows  <- as.data.frame(scores(rda_landscape, display = "bp",      scaling = 2))
# data_rda, not data_z. The two coincide now that no unit is excluded, but the
# ordination object remains the authority on which rows it used, and reading the
# genus from anywhere else is how the earlier mismatch ("replacement has 65 rows,
# data has 63") arose.
rda_site_scores$Genus <- data_rda$GENUS
stopifnot(nrow(rda_site_scores) == nrow(data_rda))

# Disambiguate biplot labels by LULC class.
# The previous compact labels collapsed forest and agropastoral predictors
# (both PD and both PROX vectors became simply "PD" and "PROX"), which made
# the figure unreadable. Forest predictors now receive the suffix "-F",
# agropastoral predictors receive "-A", and the environmental covariates keep
# their original names. Underscores are kept only to separate the metric base
# from its units (e.g. "ENN-F_m", "AREA-F_ha").
rda_biplot_arrows$Label <- rownames(rda_biplot_arrows) %>%
  # forest predictors -> "-F"
  gsub("_Forest", "-F", .) %>%
  # agropastoral predictors -> "-A"
  gsub("_Agropecuaria", "-A", .) %>%
  # drop the "Mean_" prefix for shape, ENN and AREA
  gsub("^Mean_", "", .)

axis_lab <- function(i) {
  pct <- round(100 * rda_landscape$CCA$eig[i] /
                 sum(c(rda_landscape$CCA$eig, rda_landscape$CA$eig)), 1)
  sprintf("RDA%d (%.1f%%)", i, pct)
}

genera_present <- unique(as.character(rda_site_scores$Genus))
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
colors_subset <- colors_genus[names(colors_genus) %in% genera_present]

p_rda <- ggplot() +
  geom_hline(yintercept = 0, colour = "grey85", linewidth = 0.3) +
  geom_vline(xintercept = 0, colour = "grey85", linewidth = 0.3) +
  geom_point(data = rda_site_scores,
             aes(x = RDA1, y = RDA2, colour = Genus, shape = Genus),
             size = 2.4, alpha = 0.85) +
  geom_segment(data = rda_biplot_arrows,
               aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
               arrow = arrow(length = unit(0.20, "cm")),
               colour = "grey10", linewidth = 0.7, alpha = 0.85) +
  geom_segment(data = rda_species_scores,
               aes(x = 0, y = 0, xend = RDA1, yend = RDA2),
               arrow = arrow(length = unit(0.15, "cm"), type = "closed"),
               colour = "#D55E00", linewidth = 0.5, alpha = 0.7,
               linetype = "dashed") +
  ggrepel::geom_text_repel(data = rda_biplot_arrows,
                           aes(x = RDA1, y = RDA2, label = Label),
                           size = 3, fontface = "bold",
                           segment.colour = "grey50",
                           max.overlaps = Inf) +
  ggrepel::geom_text_repel(data = rda_species_scores %>%
                             rownames_to_column("LULC"),
                           aes(x = RDA1, y = RDA2, label = LULC),
                           size = 3, fontface = "italic", colour = "#D55E00",
                           max.overlaps = Inf) +
  scale_colour_manual(values = colors_subset) +
  scale_shape_manual(values  = setNames(seq_along(genera_present), genera_present)) +
  labs(x = axis_lab(1), y = axis_lab(2)) +
  theme_classic(base_size = 10) +
  theme(plot.title    = element_text(face = "bold", size = 10, hjust = 0),
        legend.position = "right")

ggsave(filename    = "Figure_RDA_Biplot.tiff",
       plot        = p_rda,
       path        = output_path,
       width       = 174, height = 140, units = "mm",
       dpi         = 600, device = "tiff", compression = "lzw", bg = "white")
message("   Figure_RDA_Biplot.tiff exported.")


# ==============================================================================
# 12. SESSION INFO AND OUTPUT INVENTORY                                [9/9]
# ==============================================================================
message("\n[9/9] All Part II analyses completed.")
message("\nOutputs in: ", output_path)
message("  Figures:")
message("    FigS1_Collinearity_Pearson.tiff")
message("    Figure_VP_Venn.tiff")
message("    Figure_RDA_Biplot.tiff")
message("  Tables:")
message("    TableII_1_VIF_VP.csv")
message("    TableII_1_VIF_RDA.csv")
message("    TableII_1c_Correlation_PD_FRAC.csv   (recomputed on n = 67)")
message("    TableII_0_Sample_Sizes.csv           (n per analysis; all = 67)")
message("    TableII_2_MEMs_Selected.csv          (if MEMs were retained)")
message("    TableII_3_VP_Fractions.csv           (identity-checked)")
message("    TableII_3_VP_Significance.csv")
message("    TableII_3b_VP_Fractions_Without_PD.csv")
message("    TableII_4_RDA_Anova_Marginal.csv     (reported in the manuscript)")
message("    TableII_4_RDA_Anova_Terms_Sequential.csv")
message("    TableII_4b_RDA_Transformation_Comparison.csv")
message("    TableII_4_RDA_Anova_Axes.csv")
message("    TableII_5_Response_Size_Dependence.csv")
message("\nNotes:")
message("  The RDA response is Hellinger-transformed class area. The square-root")
message("  of the absolute class area used before 22 Aug 2026 carried the area of")
message("  the sampling unit into the ordination; see Section 11 and")
message("  TableII_5_Response_Size_Dependence.csv.")
message("  Term tests reported in the manuscript are MARGINAL. The sequential")
message("  decomposition is exported beside them for description only.")
message("  VIF Strategies 2 and 3 are available as commented blocks in Section 6.")
message("  Activate by uncommenting and reassigning X1 / predictors_VP before")
message("  running Section 9 (VP) and Section 11 (RDA).")

message("\n[Session]")
print(sessionInfo())
writeLines(capture.output(sessionInfo()),
           file.path(output_path, "session_info_part02.txt"))

# ==============================================================================
# END OF PART II SCRIPT (v4)
# ==============================================================================
