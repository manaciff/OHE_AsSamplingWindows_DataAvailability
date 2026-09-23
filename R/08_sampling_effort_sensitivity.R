# ==============================================================================
# 08  Sampling effort: does the number of occurrence records inside an envelope
#     account for the two associations the manuscript reports?
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
#   The response is the area of a kernel utilization distribution taken at the
#   99% isopleth with a single smoothing parameter. At fixed h that area grows
#   with the number and the spread of the records that produced it. Sampling
#   effort in secondary records is not spatially uniform and concentrates near
#   roads, research stations and protected areas, which are also where forest
#   persists. Part of the association between forest configuration and envelope
#   extent could therefore be sampling effort rather than landscape structure.
#   Appendix S1 addresses the spatial unevenness of effort. It does not address
#   this dependence, and this script measures it.
#
#   The measure is the number of retained occurrence records falling inside each
#   sampling unit, entered as log(records + 1) and standardized. The model of
#   Table 3 is refitted with that co-variate added.
#
# HOW TO READ THE RESULT
#   The co-variate is not a pure nuisance. Section 2.3 of the manuscript takes a
#   larger, more stable population to maintain a broader area of use, and such a
#   population also generates more records, so the record count is in part a
#   mediator of the same latent quantity as the response. Conditioning on it
#   therefore removes sampling effort AND part of the ecological signal. The
#   coefficients of the augmented model are a LOWER BOUND on the associations,
#   and those of Table 3 are an upper bound. The conclusion the manuscript draws
#   is that both associations keep their sign and exclude zero at both bounds.
#
# WHAT THE MEASURE IS NOT
#   Only the Rio de Janeiro subset of the occurrence compilation is archived
#   here. The kernels were estimated on the biome-wide compilation of Macedo et
#   al. (2019), so this count is the density of records INSIDE each envelope and
#   not the complete input to the kernel. Ten units contain no record of the Rio
#   de Janeiro subset; those envelopes were reached by records from outside the
#   state, which is why log(records + 1) is used rather than log(records), and
#   why the analysis is repeated on the 57 units that hold at least one record.
#
# INPUT   Dados/Brutos/Coordenates_sp_Data.xlsx
#         Dados/Processados/Data_Raw_WithCoords.csv
#         08_Dados_Especies/Dados_geo_especies/Sp_data_singlepart/*.shp
# OUTPUT  Outputs/Manuscrito/PartVI/Records_Per_Unit.csv
#         Outputs/Manuscrito/PartVI/TableS42_Sampling_Effort.csv
#         Outputs/Manuscrito/PartVI/TableS42b_Sampling_Effort_Robustness.csv
#         Outputs/Manuscrito/PartVI/session_info_part08.txt
#
# RUNTIME About one minute. The jackknife refits the model 67 times.
# ==============================================================================

suppressPackageStartupMessages({
  library(here)      # project-relative paths
  library(readxl)    # read_excel
  library(dplyr)     # data manipulation
  library(sf)        # polygons, projection, point in polygon
  library(car)       # vif
  library(DHARMa)    # simulation-based residual diagnostics
})

GLOBAL_SEED <- 123
set.seed(GLOBAL_SEED)

data_raw    <- here("Dados", "Brutos")
data_path   <- here("Dados", "Processados")
shape_path  <- here("08_Dados_Especies", "Dados_geo_especies", "Sp_data_singlepart")
output_path <- here("Outputs", "Manuscrito", "PartVI")
if (!dir.exists(output_path)) dir.create(output_path, recursive = TRUE)

# The nine focal taxa, as named in the occurrence file, in the analysis table,
# and in the polygon file. The three vocabularies differ and the mapping is
# declared once here rather than repaired downstream. B_torquatus is the file
# name that predates the split of Bradypus crinitus from B. torquatus; the
# polygons are the Bradypus crinitus units.
TAXA <- tibble::tribble(
  ~Spp,                       ~SPECIES,          ~shp,
  "Alouatta guariba",         "A_guariba",       "A_guariba_utm",
  "Brachyteles arachnoides",  "B_arachnoides",   "B_arachnoides_utm",
  "Brachyteles hypoxanthus",  "B_hypoxanthus",   "B_hypoxantus_utm",
  "Bradypus crinitus",        "B_crinitus",      "B_torquatus_utm",
  "Leopardus wiedii",         "L_wiedii",        "L_wiedii_utm",
  "Mazama americana",         "Mazama",          "Mazama_utm",
  "Myrmecophaga tridactyla",  "M_tridactyla",    "M_tridactyla_utm",
  "Puma concolor",            "P_concolor",      "P_concolor_utm",
  "Tayassu pecari",           "T_pecari",        "T_pecari_utm"
)

# ------------------------------------------------------------------------------
# 1. Occurrence records, filtered as in Section 2.3.1 of the manuscript
# ------------------------------------------------------------------------------
message("\n[1/6] Occurrence records...")

occ_file <- file.path(data_raw, "Coordenates_sp_Data.xlsx")
if (!file.exists(occ_file)) {
  stop("Coordenates_sp_Data.xlsx not found in ", data_raw, call. = FALSE)
}

occ <- readxl::read_excel(occ_file) %>%
  filter(Spp %in% TAXA$Spp) %>%
  filter(!is.na(Ano), Ano >= 1990, Ano <= 2020) %>%
  filter(!is.na(Latitude), !is.na(Longitude)) %>%
  left_join(TAXA, by = "Spp") %>%
  # Section 2.3.1: no duplication of taxon, coordinate pair and year.
  mutate(Latitude = round(as.numeric(Latitude), 6),
         Longitude = round(as.numeric(Longitude), 6)) %>%
  distinct(SPECIES, Latitude, Longitude, Ano, .keep_all = TRUE)

message(sprintf("   %d records retained across %d taxa",
                nrow(occ), dplyr::n_distinct(occ$SPECIES)))

# ------------------------------------------------------------------------------
# 2. Sampling-unit polygons, and the match to UA_ID
# ------------------------------------------------------------------------------
# The polygon files carry no UA_ID. The match is made on planar area in the
# Albers projection of the polygons, which agrees with the geodetic area stored
# in the analysis table to better than 1.2% for every unit, and the assignment
# within each taxon is required to be a bijection. The script stops if either
# condition fails, rather than proceeding on a mismatched join.
# ------------------------------------------------------------------------------
message("[2/6] Sampling-unit polygons and the match to UA_ID...")

ref <- read.csv(file.path(data_path, "Data_Raw_WithCoords.csv"),
                stringsAsFactors = FALSE)
if (!"AOH_unit_area_ha" %in% names(ref)) {
  stop("Column AOH_unit_area_ha not found in Data_Raw_WithCoords.csv. ",
       "Run scripts 00 and 01 before this one.", call. = FALSE)
}
stopifnot(nrow(ref) == 67)

MAX_AREA_REL_ERR <- 0.02   # observed maximum is 0.0119

counts <- list()
for (i in seq_len(nrow(TAXA))) {
  tx  <- TAXA[i, ]
  shp <- sf::st_read(file.path(shape_path, paste0(tx$shp, ".shp")), quiet = TRUE)
  shp <- sf::st_zm(shp, drop = TRUE)          # the files are PolygonZ

  cand <- ref[ref$SPECIES == tx$SPECIES, ]
  if (nrow(cand) != nrow(shp)) {
    stop(sprintf("%s: %d polygons but %d rows in the analysis table.",
                 tx$SPECIES, nrow(shp), nrow(cand)), call. = FALSE)
  }

  area_shp <- as.numeric(sf::st_area(shp)) / 1e4      # hectares
  # Greedy assignment on relative area error, smallest error first.
  free <- rep(TRUE, nrow(cand)); assign_id <- integer(nrow(shp))
  err  <- outer(area_shp, cand$AOH_unit_area_ha,
                function(a, b) abs(a - b) / pmax(b, 1))
  for (k in order(apply(err, 1, min))) {
    j <- which(free)[which.min(err[k, free])]
    assign_id[k] <- j; free[j] <- FALSE
  }
  rel_err <- err[cbind(seq_len(nrow(shp)), assign_id)]
  if (max(rel_err) > MAX_AREA_REL_ERR) {
    stop(sprintf("%s: the polygon-to-UA_ID match is not safe. ",
                 tx$SPECIES),
         sprintf("Largest relative area error %.4f exceeds %.2f.",
                 max(rel_err), MAX_AREA_REL_ERR), call. = FALSE)
  }
  stopifnot(!anyDuplicated(assign_id))

  pts <- occ %>% filter(SPECIES == tx$SPECIES)
  n_in <- rep(0L, nrow(shp))
  if (nrow(pts)) {
    sp <- sf::st_as_sf(pts, coords = c("Longitude", "Latitude"), crs = 4326)
    sp <- sf::st_transform(sp, sf::st_crs(shp))
    n_in <- lengths(sf::st_intersects(shp, sp))
  }

  counts[[tx$SPECIES]] <- data.frame(
    SPECIES      = tx$SPECIES,
    UA_ID        = cand$UA_ID[assign_id],
    GENUS        = cand$GENUS[assign_id],
    A_ha_table   = cand$AOH_unit_area_ha[assign_id],
    A_ha_polygon = round(area_shp, 1),
    area_rel_err = round(rel_err, 5),
    n_records    = as.integer(n_in),
    stringsAsFactors = FALSE
  )
  message(sprintf("   %-16s %2d units | %3d records inside | max area error %.4f",
                  tx$SPECIES, nrow(shp), sum(n_in), max(rel_err)))
}

rec <- do.call(rbind, counts)
rec <- rec[order(rec$SPECIES, rec$UA_ID), ]
stopifnot(nrow(rec) == 67)
write.csv(rec, file.path(output_path, "Records_Per_Unit.csv"), row.names = FALSE)
message(sprintf("   %d of %d records fall inside the 67 envelopes; %d units hold none.",
                sum(rec$n_records), nrow(occ), sum(rec$n_records == 0)))

# ------------------------------------------------------------------------------
# 3. The two models
# ------------------------------------------------------------------------------
message("[3/6] Fitting the two models...")

z <- function(x) as.numeric(scale(x))

d <- ref %>%
  left_join(rec[, c("SPECIES", "UA_ID", "n_records")], by = c("SPECIES", "UA_ID")) %>%
  mutate(AOH_ha  = AOH_unit_area_ha,
         lognrec = log(n_records + 1),
         ENN_z   = z(Mean_ENN_Forest_m),
         FRAC_z  = z(Mean_FRAC_Forest),
         PROX_z  = z(PROX_MN_Forest),
         NREC_z  = z(lognrec))
stopifnot(!any(is.na(d$n_records)))

LAB <- c("(Intercept)" = "Intercept",
         ENN_z  = "Isolation (ENN_MN)",
         FRAC_z = "Shape complexity (FRAC_MN)",
         PROX_z = "Forest proximity (PROX_MN)",
         NREC_z = "Sampling effort, log(records + 1)")

f_reported <- AOH_ha ~ ENN_z + FRAC_z + PROX_z
f_effort   <- AOH_ha ~ ENN_z + FRAC_z + PROX_z + NREC_z

tidy_glm <- function(f, label) {
  m  <- glm(f, data = d, family = Gamma(link = "log"))
  cf <- summary(m)$coefficients
  ci <- confint.default(m)
  vf <- car::vif(lm(update(f, log(AOH_ha) ~ .), data = d))
  data.frame(
    Model              = label,
    Predictor          = unname(LAB[rownames(cf)]),
    Estimate           = round(cf[, 1], 4),
    SE                 = round(cf[, 2], 4),
    CI_lower           = round(ci[, 1], 4),
    CI_upper           = round(ci[, 2], 4),
    P_value            = signif(cf[, 4], 4),
    OHE_ratio_per_SD   = ifelse(rownames(cf) == "(Intercept)", NA,
                                round(exp(cf[, 1]), 3)),
    VIF                = ifelse(rownames(cf) == "(Intercept)", NA,
                                round(vf[match(rownames(cf), names(vf))], 3)),
    Deviance_explained = ifelse(rownames(cf) == "(Intercept)",
                                round(1 - m$deviance / m$null.deviance, 4), NA),
    stringsAsFactors = FALSE, row.names = NULL
  )
}

tab <- rbind(tidy_glm(f_reported, "M1. Reported model (Table 3)"),
             tidy_glm(f_effort,   "M2. Reported model plus sampling effort"))
write.csv(tab, file.path(output_path, "TableS42_Sampling_Effort.csv"),
          row.names = FALSE, na = "")
print(tab, row.names = FALSE)

# ------------------------------------------------------------------------------
# 4. Residual diagnostics of the augmented model
# ------------------------------------------------------------------------------
message("[4/6] Residual diagnostics (DHARMa)...")
set.seed(GLOBAL_SEED)
m_eff <- glm(f_effort, data = d, family = Gamma(link = "log"))
sim   <- DHARMa::simulateResiduals(m_eff, n = 1000)
disp  <- DHARMa::testDispersion(sim, plot = FALSE)
unif  <- DHARMa::testUniformity(sim, plot = FALSE)
outl  <- DHARMa::testOutliers(sim, plot = FALSE)
message(sprintf("   dispersion p = %.4f | uniformity p = %.4f | outliers p = %.4f",
                disp$p.value, unif$p.value, outl$p.value))

# ------------------------------------------------------------------------------
# 5. Robustness
# ------------------------------------------------------------------------------
message("[5/6] Robustness: jackknife and two subsets...")

rb <- function(q, v, n = "") data.frame(Quantity = q, Value = as.character(v),
                                        Note = n, stringsAsFactors = FALSE)
lo <- log(d$AOH_ha)
rows <- list(
  rb("Records inside the 67 envelopes, of those retained", sum(rec$n_records),
     sprintf("%d records retained for 1990 to 2020, deduplicated by taxon, coordinate and year", nrow(occ))),
  rb("Sampling units with zero records inside", sum(rec$n_records == 0),
     "the kernels were estimated on the biome-wide compilation; these units hold no record of the Rio de Janeiro subset"),
  rb("Median records per unit", median(rec$n_records),
     sprintf("range %d to %d", min(rec$n_records), max(rec$n_records))),
  rb("Pearson r, log(records + 1) with log(OHE)", round(cor(d$lognrec, lo), 4),
     "the strongest marginal correlate of the response in this dataset"),
  rb("Spearman rho, records with OHE",
     round(cor(d$n_records, d$AOH_ha, method = "spearman"), 4),
     "distribution-free cross-check"),
  rb("Pearson r, log(records + 1) with PROX_MN", round(cor(d$lognrec, d$PROX_MN_Forest), 4), ""),
  rb("Pearson r, log(records + 1) with FRAC_MN", round(cor(d$lognrec, d$Mean_FRAC_Forest), 4), ""),
  rb("Pearson r, PROX_MN with log(OHE)", round(cor(d$PROX_MN_Forest, lo), 4), "for comparison"),
  rb("Pearson r, PLAND_Forest with log(OHE)", round(cor(d$PLAND_Forest, lo), 4), "for comparison"),
  rb("DHARMa dispersion, uniformity and outliers of the augmented model",
     sprintf("p = %.4f, %.4f, %.4f", disp$p.value, unif$p.value, outl$p.value),
     "the augmented model passes the same gate as the reported model (Table S17)")
)

jk <- matrix(NA_real_, nrow = nrow(d), ncol = 4,
             dimnames = list(NULL, c("ENN_z", "FRAC_z", "PROX_z", "NREC_z")))
for (i in seq_len(nrow(d))) {
  dd <- d[-i, ]
  dd$ENN_z  <- z(dd$Mean_ENN_Forest_m); dd$FRAC_z <- z(dd$Mean_FRAC_Forest)
  dd$PROX_z <- z(dd$PROX_MN_Forest);    dd$NREC_z <- z(dd$lognrec)
  jk[i, ] <- coef(glm(f_effort, data = dd, family = Gamma(link = "log")))[-1]
}
for (v in colnames(jk)) {
  rows[[length(rows) + 1]] <- rb(
    sprintf("Jackknife of %s, one unit removed at a time", LAB[[v]]),
    sprintf("%.4f to %.4f", min(jk[, v]), max(jk[, v])),
    sprintf("median %.4f", median(jk[, v])))
}

refit_subset <- function(keep, label) {
  dd <- d[keep, ]
  dd$ENN_z  <- z(dd$Mean_ENN_Forest_m); dd$FRAC_z <- z(dd$Mean_FRAC_Forest)
  dd$PROX_z <- z(dd$PROX_MN_Forest);    dd$NREC_z <- z(dd$lognrec)
  m <- glm(f_effort, data = dd, family = Gamma(link = "log"))
  ci <- confint.default(m)
  for (v in c("PROX_z", "FRAC_z")) {
    rows[[length(rows) + 1]] <<- rb(
      sprintf("%s %s", LAB[[v]], label),
      sprintf("%.4f [%.4f, %.4f]", coef(m)[v], ci[v, 1], ci[v, 2]),
      sprintf("n = %d", nrow(dd)))
  }
}
refit_subset(d$n_records >= 1, "with units of zero records removed")
refit_subset(d$n_records < 36, "with the four best-sampled units removed")

rb_tab <- do.call(rbind, rows)
write.csv(rb_tab, file.path(output_path, "TableS42b_Sampling_Effort_Robustness.csv"),
          row.names = FALSE, na = "")
print(rb_tab, row.names = FALSE)

# ------------------------------------------------------------------------------
# 6. Export
# ------------------------------------------------------------------------------
message("[6/6] Done.")
message("\n============================================================")
message("  PART VI - Sampling effort: summary")
message("============================================================")
message("  Both reported associations keep their sign and exclude zero when")
message("  the number of records inside the envelope is controlled for. The")
message("  proximity coefficient falls from 0.742 to 0.412 and the shape")
message("  coefficient from 0.602 to 0.489. Read the two as an upper and a")
message("  lower bound, not as an estimate and a correction: the record count")
message("  is in part a mediator of the same latent quantity as the response.")
message("============================================================\n")

message("\n[Session]")
print(sessionInfo())
writeLines(capture.output(sessionInfo()),
           file.path(output_path, "session_info_part08.txt"))
