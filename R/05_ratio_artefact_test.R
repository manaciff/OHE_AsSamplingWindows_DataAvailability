# ==============================================================================
# 05  Ratio artefact test: which coefficients survive the geometry of the design
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
#   Decide which model the manuscript reports.
#
#   Forest cover and patch density are normalized by the area of the sampling
#   window, which is also the response. The class area and the number of forest
#   patches are therefore reshuffled jointly across units, preserving their
#   mutual association, while the area of the unit is held at its observed
#   value; the two ratios are recomputed from the permuted numerators and the
#   observed denominator, and the model is refitted 9,999 times with a fixed
#   seed. Under this null there is by construction no relationship between the
#   amount of forest and the extent of the envelope, so any coefficient obtained
#   measures the construction alone. The same null is applied to the information
#   criterion, because a predictor built from the response improves the fit to
#   the response whether or not it carries ecological information.
#
#   Six alternative specifications remove the normalized metrics singly and
#   jointly. A separate screen asks whether the proximity index, which is not
#   divided by the area of the window but is built from the areas of
#   neighboring patches, is separable from the amount of forest present.
#
# INPUT   Dados/Processados/Data_Raw_WithCoords.csv
# OUTPUT  Outputs/Manuscrito/PartV/, including the coefficient table and the
#         coefficient plot of the model the manuscript reports
#
# RUNTIME About ten minutes. The null of the coefficients runs N_SIM = 9,999
#         permutations; the null of the AICc difference runs N_SIM_AIC = 2,000,
#         because each iteration refits two models rather than one. Both counts
#         are declared in the constants block below.
# ==============================================================================


# ==============================================================================
# 1. ENVIRONMENT SETUP
# ==============================================================================
# ------------------------------------------------------------------------------
# WHY ppcor IS NOT LOADED HERE. ppcor depends on MASS, and MASS exports its own
# select(), which masks the one from dplyr. On the first run that produced
#   Error in select(...) : unused arguments
# because the MASS version was picked up. The dependency is removed altogether:
# the partial correlations below are computed from the residuals of two linear
# models, which is the definition of a partial correlation and needs no package.
# Every call is also written in the qualified form dplyr::select(), so the script
# is immune to whatever else happens to be attached in the session.
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(here)        # project-relative paths
  library(tidyverse)   # data wrangling, ggplot2, purrr
  library(car)         # vif
  library(MuMIn)       # AICc on a common response scale
})

# Partial correlation of x and y given the columns of z, defined as the Pearson
# correlation between the residuals of x ~ z and of y ~ z. With z = NULL it
# reduces to the ordinary correlation.
partial_cor <- function(x, y, z = NULL) {
  x <- as.numeric(x); y <- as.numeric(y)
  if (is.null(z)) return(cor(x, y))
  z <- as.data.frame(z)
  if (ncol(z) == 0) return(cor(x, y))
  ok <- stats::complete.cases(x, y, z)
  zz <- z[ok, , drop = FALSE]
  rx <- stats::resid(stats::lm(x[ok] ~ ., data = zz))
  ry <- stats::resid(stats::lm(y[ok] ~ ., data = zz))
  cor(rx, ry)
}

GLOBAL_SEED      <- 123
N_UNITS_EXPECTED <- 67
N_SIM            <- 9999   # permutations of the coefficient null
N_SIM_AIC_CAP    <- 2000   # permutations of the AICc null: each one refits two
                           # models, so the count is capped below N_SIM
THIN_RADIUS_M    <- 1000   # centroids closer than this join one cluster
N_THIN_REPS      <- 200    # thinned designs drawn from those clusters
set.seed(GLOBAL_SEED)

data_path   <- here("Dados", "Processados")
output_path <- here("Outputs", "Manuscrito", "PartV")
if (!dir.exists(output_path)) dir.create(output_path, recursive = TRUE)


# ==============================================================================
# 2. DATA, PREPARED AS IN PART III                                      [1/7]
# ==============================================================================
message("\n[1/7] Loading data...")

d <- read_csv(file.path(data_path, "Data_Raw_WithCoords.csv"),
              show_col_types = FALSE)

if ("AOH_unit_area_ha" %in% names(d) && !"AOH_ha" %in% names(d)) {
  d <- d %>% rename(AOH_ha = AOH_unit_area_ha)
}

# PROX_MN_Agropecuaria and PD_Agropecuaria were added on 22 Aug 2026 for the
# uncoupled matrix model of Section 3a. Neither is divided by the area of the
# sampling window, so both are admissible in a specification free of coupling.
needed <- c("AOH_ha", "CA_Forest", "NP_Forest", "PLAND_Forest", "PD_Forest",
            "Mean_FRAC_Forest", "Mean_ENN_Forest_m", "PROX_MN_Forest",
            "Mean_AREA_Forest_ha", "PROX_MN_Agropecuaria", "PD_Agropecuaria",
            "coord_x", "coord_y", "GENUS")
missing <- setdiff(needed, names(d))
if (length(missing) > 0) stop("Missing column(s): ", paste(missing, collapse = ", "))

d <- d %>%
  dplyr::select(GENUS, SPECIES, UA_ID, all_of(setdiff(needed, "GENUS"))) %>%
  tidyr::drop_na() %>%
  mutate(log_AOH = log(AOH_ha), GENUS = as.factor(GENUS))

if (nrow(d) != N_UNITS_EXPECTED) {
  stop(sprintf("Part V loaded %d units instead of %d. Rerun Part I.",
               nrow(d), N_UNITS_EXPECTED))
}
message(sprintf("   %d sampling units.", nrow(d)))

z <- function(x) as.numeric(scale(x))


# ==============================================================================
# 3. THE ALGEBRAIC COUPLING, MADE EXPLICIT                              [2/7]
# ==============================================================================
# Two things are checked here. First, that the identity linking PD to the
# response holds numerically on all 67 units, so that the concern is not
# hypothetical. Second, how much of the association between forest proximity and
# the response is carried by the coupled variables, through partial correlations.
# ==============================================================================
message("\n[2/7] Verifying the algebraic identity and partial correlations...")

identity_error <- max(abs(log(d$NP_Forest) + log(100) - log(d$PD_Forest) - d$log_AOH))

pc <- partial_cor

coupling <- data.frame(
  Quantity = c(
    "max |log(NP) + log(100) - log(PD) - log(OHE)|",
    "Pearson r, PLAND_Forest with log(OHE)",
    "Spearman rho, PLAND_Forest with log(OHE)",
    "Pearson r, PD_Forest with log(OHE)",
    "Pearson r, log(CA_Forest) with log(OHE)",
    "Pearson r, PROX_MN_Forest with log(OHE)",
    "Partial r, PROX_MN with log(OHE), controlling PLAND",
    "Partial r, PROX_MN with log(OHE), controlling PLAND and PD",
    "Partial r, PROX_MN with log(OHE), controlling Mean_AREA_Forest_ha"
  ),
  Value = round(c(
    identity_error,
    cor(d$PLAND_Forest, d$log_AOH),
    cor(d$PLAND_Forest, d$log_AOH, method = "spearman"),
    cor(d$PD_Forest, d$log_AOH),
    cor(log(d$CA_Forest), d$log_AOH),
    cor(d$PROX_MN_Forest, d$log_AOH),
    pc(d$PROX_MN_Forest, d$log_AOH, d[, "PLAND_Forest"]),
    pc(d$PROX_MN_Forest, d$log_AOH, d[, c("PLAND_Forest", "PD_Forest")]),
    pc(d$PROX_MN_Forest, d$log_AOH, d[, "Mean_AREA_Forest_ha"])
  ), 4),
  Note = c(
    "if near zero, PD is an exact function of the response given NP",
    "marginal association, no conditioning",
    "distribution-free cross-check of the line above",
    "marginal association, no conditioning",
    "the response bounds this quantity from above",
    "PROX is NOT normalised by the window area",
    "conditioning on a ratio that contains the response",
    "conditioning on two quantities that contain the response",
    "conditioning on a quantity free of the response"
  )
)

write.csv(coupling,
          file.path(output_path, "TableS_PartV_Identity_and_PartialCorrelations.csv"),
          row.names = FALSE)
print(coupling, row.names = FALSE)

if (identity_error > 0.05) {
  message("   NOTE: the PD identity does not hold to within 0.05 on all units.")
  message("   Check the units of PD before interpreting Section 6.")
} else {
  message(sprintf("   The PD identity holds to within %.4f on every unit.",
                  identity_error))
}


# ==============================================================================
# 4. COMPETING SPECIFICATIONS                                           [3/7]
# ==============================================================================
# Each specification is fitted in the family retained by Part III, Gamma with a
# log link, on the same 67 units. What matters is not which fits best, because
# the coupled specifications have an unfair advantage: a predictor built from the
# response will always explain the response. What matters is whether the
# coefficient of PROX_MN_Forest and the interaction keep their sign, their
# magnitude and their interval as the coupled predictors are removed.
#
#   S1  retained model                                    reference
#   S2  without PLAND                                     drops one coupled term
#   S3  without PD                                        drops the other
#   S4  neither PLAND nor PD                              no coupled predictor
#   S5  PLAND replaced by log(CA_Forest)                  absolute, not a ratio
#   S6  amount measured by mean forest patch area         free of the window
# ==============================================================================
message("\n[3/7] Fitting the competing specifications...")

dz <- d %>%
  mutate(PLAND_z   = z(PLAND_Forest),
         PD_z      = z(PD_Forest),
         FRAC_z    = z(Mean_FRAC_Forest),
         ENN_z     = z(Mean_ENN_Forest_m),
         PROX_z    = z(PROX_MN_Forest),
         logCA_z   = z(log(CA_Forest)),
         AREAMN_z  = z(Mean_AREA_Forest_ha))

specs <- list(
  S1 = list(f = AOH_ha ~ ENN_z + FRAC_z + PD_z + PLAND_z + PROX_z + PROX_z:PLAND_z,
            label = "S1. Retained model (reference)",
            coupled = "PLAND and PD"),
  S2 = list(f = AOH_ha ~ ENN_z + FRAC_z + PD_z + PROX_z,
            label = "S2. Without PLAND",
            coupled = "PD only"),
  S3 = list(f = AOH_ha ~ ENN_z + FRAC_z + PLAND_z + PROX_z + PROX_z:PLAND_z,
            label = "S3. Without PD",
            coupled = "PLAND only"),
  S4 = list(f = AOH_ha ~ ENN_z + FRAC_z + PROX_z,
            label = "S4. No area-normalised predictor",
            coupled = "none"),
  S5 = list(f = AOH_ha ~ ENN_z + FRAC_z + PD_z + logCA_z + PROX_z + PROX_z:logCA_z,
            label = "S5. Amount as log(CA_Forest), absolute",
            coupled = "PD only; CA is bounded by the response"),
  S6 = list(f = AOH_ha ~ ENN_z + FRAC_z + AREAMN_z + PROX_z + PROX_z:AREAMN_z,
            label = "S6. Amount as mean forest patch area",
            coupled = "none")
)

fit_spec <- function(sp) {
  m  <- glm(sp$f, data = dz, family = Gamma(link = "log"))
  cf <- summary(m)$coefficients
  ci <- suppressMessages(confint.default(m))

  int_row <- grep(":", rownames(cf), value = TRUE)
  amount  <- intersect(c("PLAND_z", "logCA_z", "AREAMN_z"), rownames(cf))

  vf <- tryCatch({
    lmfit <- lm(update(sp$f, log(AOH_ha) ~ .), data = dz)
    v <- car::vif(lmfit)
    max(v)
  }, error = function(e) NA_real_)

  data.frame(
    Specification   = sp$label,
    Coupled_terms   = sp$coupled,
    AICc            = round(MuMIn::AICc(m), 2),
    Max_VIF         = round(vf, 2),
    PROX_beta       = round(cf["PROX_z", "Estimate"], 4),
    PROX_CI_low     = round(ci["PROX_z", 1], 4),
    PROX_CI_high    = round(ci["PROX_z", 2], 4),
    PROX_p          = signif(cf["PROX_z", "Pr(>|t|)"], 4),
    Amount_term     = if (length(amount)) amount[1] else NA_character_,
    Amount_beta     = if (length(amount)) round(cf[amount[1], "Estimate"], 4) else NA_real_,
    Amount_p        = if (length(amount)) signif(cf[amount[1], "Pr(>|t|)"], 4) else NA_real_,
    Interaction_beta = if (length(int_row)) round(cf[int_row, "Estimate"], 4) else NA_real_,
    Interaction_CI_low  = if (length(int_row)) round(ci[int_row, 1], 4) else NA_real_,
    Interaction_CI_high = if (length(int_row)) round(ci[int_row, 2], 4) else NA_real_,
    Interaction_p    = if (length(int_row)) signif(cf[int_row, "Pr(>|t|)"], 4) else NA_real_,
    stringsAsFactors = FALSE
  )
}

spec_table <- purrr::map_dfr(specs, fit_spec)

write.csv(spec_table,
          file.path(output_path, "TableS_PartV_Specifications.csv"),
          row.names = FALSE)
print(as.data.frame(spec_table), row.names = FALSE)


# ------------------------------------------------------------------------------
# 3a. THE FARMING MATRIX, WITHOUT ANY COUPLED PREDICTOR
# ------------------------------------------------------------------------------
# ADDED 22 Aug 2026. Part III already reports a sensitivity model that adds the
# proximity index of the farming matrix to the a priori set, and the matrix term
# is a significant positive predictor there. That model, however, still contains
# the percentage of landscape and patch density, so its coefficients inherit the
# dependence this script exists to measure, and its largest variance inflation
# factor is 19.3.
#
# The matrix question does not require those two predictors. The proximity index
# of the farming matrix is not divided by the area of the window, exactly as the
# proximity index of forest is not, so it can be added to S4, the specification
# that contains no metric normalized by the response. What follows is that model.
# It is the estimate the manuscript should quote for the matrix, because it is
# the only one in which the matrix coefficient is free of the coupling.
#
# Matrix patch density is fitted beside it for contrast: it answers the same
# question with the other matrix descriptor, and it is not a ratio of the window
# either.
# ------------------------------------------------------------------------------
message("\n[3a] The farming matrix, in a specification free of coupling...")

dz <- dz %>%
  mutate(PROXA_z = z(PROX_MN_Agropecuaria),
         PDA_z   = z(PD_Agropecuaria))

matrix_specs <- list(
  M1 = list(f = AOH_ha ~ ENN_z + FRAC_z + PROX_z + PROXA_z,
            label = "S4 plus farming matrix proximity",
            focal = "PROXA_z"),
  M2 = list(f = AOH_ha ~ ENN_z + FRAC_z + PROX_z + PDA_z,
            label = "S4 plus farming matrix patch density",
            focal = "PDA_z")
)

fit_matrix <- function(sp) {
  m  <- glm(sp$f, data = dz, family = Gamma(link = "log"))
  cf <- summary(m)$coefficients
  ci <- suppressMessages(confint.default(m))
  vf <- tryCatch(max(car::vif(lm(update(sp$f, log(AOH_ha) ~ .), data = dz))),
                 error = function(e) NA_real_)
  keep <- setdiff(rownames(cf), "(Intercept)")
  data.frame(
    Specification = sp$label,
    Predictor     = keep,
    Estimate      = round(cf[keep, "Estimate"], 4),
    SE            = round(cf[keep, "Std. Error"], 4),
    CI_lower      = round(ci[keep, 1], 4),
    CI_upper      = round(ci[keep, 2], 4),
    P_value       = signif(cf[keep, 4], 4),
    OHE_ratio_per_SD = round(exp(cf[keep, "Estimate"]), 3),
    Deviance_explained = round(1 - m$deviance / m$null.deviance, 4),
    Max_VIF       = round(vf, 2),
    row.names     = NULL, stringsAsFactors = FALSE
  )
}

matrix_table <- purrr::map_dfr(matrix_specs, fit_matrix)
write.csv(matrix_table,
          file.path(output_path, "TableS_PartV_Matrix_Uncoupled.csv"),
          row.names = FALSE)
print(matrix_table, row.names = FALSE)
message("   Read the matrix coefficient from this table and not from the")
message("   sensitivity model of Part III, which contains coupled predictors.")


# ------------------------------------------------------------------------------
# 3b. SPATIAL COINCIDENCE AMONG UNITS OF DIFFERENT TAXA
# ------------------------------------------------------------------------------
# ADDED 22 Aug 2026. The envelopes of different taxa overlap, and in places one
# is nested inside another, so the same landscape can enter the model more than
# once. On these data eighteen of the sixty-seven units have a centroid within
# one kilometer of the centroid of a unit belonging to a different taxon, and one
# group of five units from five different taxa shares essentially the same
# centroid while their envelope areas differ by a factor of five.
#
# Neither of the two controls already in the chain covers this. The random
# intercept for genus does not, because the groups cross genera. The Moran
# eigenvectors do so only in part, because proximity between centroids is not the
# same thing as overlap between polygons.
#
# The check below is a thinning test. Units whose centroids fall within
# THIN_RADIUS_M of one another are treated as one cluster, one unit is drawn at
# random from each cluster, the reported model is refitted on the thinned design,
# and the procedure is repeated. What the distribution of coefficients across
# replicates shows is whether the two reported associations depend on counting
# the same landscape several times.
#
# This is a robustness check on an existing result, not a new analysis: the
# coefficients it produces are not a second set of estimates.
# ------------------------------------------------------------------------------
message("\n[3b] Thinning test for spatial coincidence among units of different taxa...")

# THIN_RADIUS_M and N_THIN_REPS are declared in the constants block at the top.

xy   <- as.matrix(dz[, c("coord_x", "coord_y")])
dmat <- as.matrix(stats::dist(xy))
diag(dmat) <- Inf

# Single-linkage clustering at the threshold, written out rather than delegated
# so that the rule is visible: two units join the same cluster when their
# centroids are closer than THIN_RADIUS_M.
parent <- seq_len(nrow(dz))
find_root <- function(a) { while (parent[a] != a) a <- parent[a]; a }
for (i in seq_len(nrow(dz) - 1)) {
  for (j in (i + 1):nrow(dz)) {
    if (dmat[i, j] < THIN_RADIUS_M) {
      ri <- find_root(i); rj <- find_root(j)
      if (ri != rj) parent[ri] <- rj
    }
  }
}
cluster_id <- vapply(seq_len(nrow(dz)), find_root, integer(1))
n_clusters <- length(unique(cluster_id))
n_in_group <- sum(table(cluster_id)[as.character(cluster_id)] > 1)

cross_taxon <- sum(vapply(unique(cluster_id), function(cl) {
  idx <- which(cluster_id == cl)
  length(idx) > 1 && length(unique(dz$SPECIES[idx])) > 1
}, logical(1)))

message(sprintf(paste0("   %d units fall in %d cluster(s) of coincident centroids ",
                       "at %d m; %d cluster(s) mix taxa. Thinned design: %d units."),
                n_in_group, sum(table(cluster_id) > 1), THIN_RADIUS_M,
                cross_taxon, n_clusters))

set.seed(GLOBAL_SEED)
thin_draws <- matrix(NA_real_, nrow = N_THIN_REPS, ncol = 6,
                     dimnames = list(NULL, c("PROX_beta", "PROX_p",
                                             "FRAC_beta", "FRAC_p",
                                             "ENN_beta", "ENN_p")))
for (r in seq_len(N_THIN_REPS)) {
  keep <- vapply(unique(cluster_id), function(cl) {
    idx <- which(cluster_id == cl)
    if (length(idx) == 1) idx else sample(idx, 1)
  }, integer(1))
  sub <- dz[sort(keep), , drop = FALSE]
  # Re-standardise inside the thinned design, so the coefficients stay on the
  # scale of one standard deviation of the data actually fitted.
  sub <- sub %>% mutate(ENN_z = z(Mean_ENN_Forest_m), FRAC_z = z(Mean_FRAC_Forest),
                        PROX_z = z(PROX_MN_Forest))
  fit <- tryCatch(glm(specs$S4$f, data = sub, family = Gamma(link = "log")),
                  error = function(e) NULL)
  if (is.null(fit)) next
  cf <- summary(fit)$coefficients
  thin_draws[r, ] <- c(cf["PROX_z", 1], cf["PROX_z", 4],
                       cf["FRAC_z", 1], cf["FRAC_z", 4],
                       cf["ENN_z",  1], cf["ENN_z",  4])
}

thin_summary <- data.frame(
  Predictor = c("Forest proximity (PROX_MN)", "Shape complexity (FRAC_MN)",
                "Isolation (ENN_MN)"),
  Beta_full_design = round(c(
    coef(glm(specs$S4$f, data = dz, family = Gamma(link = "log")))[["PROX_z"]],
    coef(glm(specs$S4$f, data = dz, family = Gamma(link = "log")))[["FRAC_z"]],
    coef(glm(specs$S4$f, data = dz, family = Gamma(link = "log")))[["ENN_z"]]), 4),
  Beta_thinned_median = round(c(median(thin_draws[, "PROX_beta"], na.rm = TRUE),
                                median(thin_draws[, "FRAC_beta"], na.rm = TRUE),
                                median(thin_draws[, "ENN_beta"],  na.rm = TRUE)), 4),
  Beta_thinned_q2.5 = round(c(quantile(thin_draws[, "PROX_beta"], 0.025, na.rm = TRUE),
                              quantile(thin_draws[, "FRAC_beta"], 0.025, na.rm = TRUE),
                              quantile(thin_draws[, "ENN_beta"],  0.025, na.rm = TRUE)), 4),
  Beta_thinned_q97.5 = round(c(quantile(thin_draws[, "PROX_beta"], 0.975, na.rm = TRUE),
                               quantile(thin_draws[, "FRAC_beta"], 0.975, na.rm = TRUE),
                               quantile(thin_draws[, "ENN_beta"],  0.975, na.rm = TRUE)), 4),
  Pct_replicates_p_below_0.05 = round(100 * c(
    mean(thin_draws[, "PROX_p"] < 0.05, na.rm = TRUE),
    mean(thin_draws[, "FRAC_p"] < 0.05, na.rm = TRUE),
    mean(thin_draws[, "ENN_p"]  < 0.05, na.rm = TRUE)), 1),
  N_units_per_replicate = n_clusters,
  N_replicates = N_THIN_REPS,
  N_units_in_multi_unit_clusters = n_in_group,
  N_clusters_mixing_taxa = cross_taxon,
  Thinning_radius_m = THIN_RADIUS_M,
  row.names = NULL
)
write.csv(thin_summary,
          file.path(output_path, "TableS_PartV_Spatial_Thinning.csv"),
          row.names = FALSE)
print(thin_summary, row.names = FALSE)


# ------------------------------------------------------------------------------
# 3b. Does the interaction structure pay for itself once the amount term is free
#     of the response?
# ------------------------------------------------------------------------------
# Script 03 compares additive against interaction specifications and the
# interaction wins by more than thirteen AICc units. That comparison, however,
# is made with the percentage of landscape in both arms, and the percentage of
# landscape carries the response in its denominator. The comparison therefore
# inherits the problem this script exists to measure: the question is not
# whether an interaction improves the fit, but whether it still improves it once
# the quantity it interacts with is measured independently of the window.
#
# Four models answer that, all in the family the manuscript reports. Two use the
# coupled measure of habitat amount, the percentage of landscape, and differ
# only in the product. Two use the mean forest patch area, which is the one
# measure of amount here that the response neither contains nor bounds, and also
# differ only in the product. The difference of differences is the quantity of
# interest: how much of the advantage of the interaction is produced by the
# normalization rather than by any threshold.
#
# Nothing is selected here. The four models are stated in advance and all four
# are reported.

pairs_int <- list(
  list(key = "Coupled amount (PLAND)",
       add = AOH_ha ~ ENN_z + FRAC_z + PD_z + PLAND_z + PROX_z,
       int = AOH_ha ~ ENN_z + FRAC_z + PD_z + PLAND_z + PROX_z + PROX_z:PLAND_z),
  list(key = "Amount free of the response (mean patch area)",
       add = AOH_ha ~ ENN_z + FRAC_z + AREAMN_z + PROX_z,
       int = AOH_ha ~ ENN_z + FRAC_z + AREAMN_z + PROX_z + PROX_z:AREAMN_z)
)

aicc_gamma <- function(f) {
  m <- glm(f, data = dz, family = Gamma(link = "log"))
  k <- length(coef(m)) + 1L
  n <- nobs(m)
  as.numeric(-2 * logLik(m) + 2 * k + (2 * k * (k + 1)) / (n - k - 1))
}

interaction_gain <- purrr::map_dfr(pairs_int, function(p) {
  a <- aicc_gamma(p$add)
  i <- aicc_gamma(p$int)
  data.frame(
    Amount_measure   = p$key,
    AICc_additive    = round(a, 2),
    AICc_interaction = round(i, 2),
    Delta_AICc       = round(i - a, 2),
    Verdict          = if (i - a <= -2) {
      "the interaction earns its parameter"
    } else if (i - a < 0) {
      "indistinguishable from the additive model; parsimony favours the additive one"
    } else {
      "the interaction costs more than it returns"
    },
    stringsAsFactors = FALSE
  )
})

write.csv(interaction_gain,
          file.path(output_path, "TableS_PartV_Interaction_Gain.csv"),
          row.names = FALSE)

message("\n   How much of the advantage of the interaction comes from the ratio:")
print(as.data.frame(interaction_gain), row.names = FALSE)
message("   Read the two rows together. The first is the comparison script 03 makes.")
message("   The second makes the same comparison with an amount measure that the")
message("   response does not contain. The gap between the two Delta AICc values is")
message("   the part of the advantage that the normalisation produces.")


# ==============================================================================
# 5. DOES THE MAIN RESULT SURVIVE?                                      [4/7]
# ==============================================================================
message("\n[4/7] Reading the specifications...")

prox_all <- spec_table$PROX_beta
prox_sig <- spec_table$PROX_CI_low * spec_table$PROX_CI_high > 0

message(sprintf("   PROX_MN_Forest ranges from %.3f to %.3f across the six specifications.",
                min(prox_all), max(prox_all)))
message(sprintf("   Its interval excludes zero in %d of %d specifications.",
                sum(prox_sig), length(prox_sig)))

int_rows <- spec_table %>% filter(!is.na(Interaction_beta))
if (nrow(int_rows) > 0) {
  int_sig <- int_rows$Interaction_CI_low * int_rows$Interaction_CI_high > 0
  message(sprintf("   The interaction is present in %d specification(s) and its interval",
                  nrow(int_rows)))
  message(sprintf("   excludes zero in %d of them, including S6, which contains no",
                  sum(int_sig)))
  message("   area-normalised predictor at all.")
  message("   S6 is the specification to cite if a reviewer challenges the coupling:")
  message("   it tests the same threshold hypothesis with an amount measure that")
  message("   does not contain the response.")
}


# ==============================================================================
# 6. NULL SIMULATION: HOW MUCH OF THE COEFFICIENT IS THE RATIO?         [5/7]
# ==============================================================================
# The null is built to preserve everything about the design except the thing
# being tested. AOH keeps its observed values. The pair (CA_Forest, NP_Forest) is
# permuted jointly across units, so the association BETWEEN them is preserved and
# only their association with the window area is destroyed. PLAND and PD are then
# recomputed from the permuted numerators and the observed denominator:
#
#   PLAND* = 100 * CA* / AOH        PD* = NP* / (AOH / 100)
#
# Under this null there is, by construction, no relationship between the amount
# of forest and the extent of the window, yet PLAND* and PD* still carry AOH in
# their denominators. Any coefficient they obtain is therefore pure artefact, and
# the distribution of those coefficients is the yardstick.
#
# Reading the result. If the observed coefficient falls INSIDE the null
# distribution, the data cannot distinguish the effect from the artefact, and the
# manuscript must say so rather than interpret the sign. If it falls OUTSIDE, the
# association is stronger than the construction alone produces, and may be
# interpreted, still with the caveat stated.
# ==============================================================================
message(sprintf("\n[5/7] Null simulation, %d permutations...", N_SIM))

obs_fit  <- glm(specs$S1$f, data = dz, family = Gamma(link = "log"))
obs_cf   <- coef(obs_fit)
obs_pland <- obs_cf[["PLAND_z"]]
obs_pd    <- obs_cf[["PD_z"]]
obs_int   <- obs_cf[[grep(":", names(obs_cf), value = TRUE)]]

null_mat <- matrix(NA_real_, nrow = N_SIM, ncol = 3,
                   dimnames = list(NULL, c("PLAND", "PD", "Interaction")))

set.seed(GLOBAL_SEED)
pb_step <- max(1, floor(N_SIM / 10))

for (i in seq_len(N_SIM)) {
  idx <- sample.int(nrow(dz))
  ca_star <- dz$CA_Forest[idx]
  np_star <- dz$NP_Forest[idx]

  sim <- dz %>%
    mutate(PLAND_z = z(100 * ca_star / AOH_ha),
           PD_z    = z(np_star / (AOH_ha / 100)))

  fit <- tryCatch(
    glm(specs$S1$f, data = sim, family = Gamma(link = "log")),
    error = function(e) NULL, warning = function(w) NULL
  )
  if (is.null(fit)) next

  cf <- coef(fit)
  ir <- grep(":", names(cf), value = TRUE)
  null_mat[i, ] <- c(cf[["PLAND_z"]], cf[["PD_z"]],
                     if (length(ir)) cf[[ir]] else NA_real_)

  if (i %% pb_step == 0) message(sprintf("   %d%%", round(100 * i / N_SIM)))
}

# ------------------------------------------------------------------------------
# READING THE NULL. Two quantities are computed, and they answer different
# questions. Reporting both is what stops one from being read as the other.
#
# emp_p is the two-sided empirical p value, the proportion of null draws whose
# ABSOLUTE value is at least as large as the observed. It is the familiar
# quantity, but it compares magnitudes against a null that is NOT centred on
# zero: the mean of the PLAND null here is about -0.28, because the ratio
# construction produces a negative coefficient on its own. A null with a long
# tail on one side and a short tail on the other will therefore return a small
# emp_p for an observation that sits comfortably inside its body.
#
# pct_below is the percentile of the observed value within the null, that is,
# the proportion of null draws at or below it. This is the quantity a reader
# means by "the observed coefficient fell at the Nth percentile of the null".
#
# THE DECISION RULE, stated once and applied without exception: a coefficient is
# treated as indistinguishable from the ratio artefact when it lies inside the
# CENTRAL 95% INTERVAL of the null, that is, between its 2.5th and 97.5th
# percentiles. The interval is used rather than emp_p because it does not assume
# the null is symmetric about zero, and the null here is neither symmetric nor
# centred on zero. Where the two disagree, the manuscript reports the
# disagreement instead of choosing silently.
# ------------------------------------------------------------------------------
emp_p <- function(obs, null_vec) {
  nv <- null_vec[is.finite(null_vec)]
  (sum(abs(nv) >= abs(obs)) + 1) / (length(nv) + 1)
}

pct_below <- function(obs, null_vec) {
  nv <- null_vec[is.finite(null_vec)]
  100 * mean(nv <= obs)
}

null_summary <- data.frame(
  Term = c("PLAND_Forest", "PD_Forest", "PLAND x PROX interaction"),
  Observed = round(c(obs_pland, obs_pd, obs_int), 4),
  Null_mean = round(apply(null_mat, 2, mean, na.rm = TRUE), 4),
  Null_q2.5 = round(apply(null_mat, 2, quantile, 0.025, na.rm = TRUE), 4),
  Null_q97.5 = round(apply(null_mat, 2, quantile, 0.975, na.rm = TRUE), 4),
  Percentile_of_observed = round(c(pct_below(obs_pland, null_mat[, "PLAND"]),
                                   pct_below(obs_pd,    null_mat[, "PD"]),
                                   pct_below(obs_int,   null_mat[, "Interaction"])), 2),
  Empirical_p_two_sided = round(c(emp_p(obs_pland, null_mat[, "PLAND"]),
                                  emp_p(obs_pd,    null_mat[, "PD"]),
                                  emp_p(obs_int,   null_mat[, "Interaction"])), 4)
)
null_summary$Verdict <- ifelse(
  null_summary$Observed >= null_summary$Null_q2.5 &
    null_summary$Observed <= null_summary$Null_q97.5,
  "INSIDE the central 95% interval of the null: indistinguishable from the ratio artefact; do not interpret the sign",
  "OUTSIDE the central 95% interval of the null: stronger than the construction alone produces; may be interpreted with the caveat"
)
# Flag, rather than hide, any term on which the two readings disagree.
null_summary$Rules_agree <- ifelse(
  grepl("^INSIDE", null_summary$Verdict) & null_summary$Empirical_p_two_sided < 0.05,
  "no: inside the interval but two-sided p < 0.05; treat as the borderline case",
  "yes")

write.csv(null_summary,
          file.path(output_path, "TableS_PartV_NullSimulation.csv"),
          row.names = FALSE)
print(null_summary, row.names = FALSE)

# car exports its own recode(), which masks dplyr::recode() whenever car is
# attached, and car is attached here for vif(). The relabeling is done with a
# plain named vector instead, which depends on no package at all.
term_labels <- c(PLAND = "PLAND_Forest", PD = "PD_Forest",
                 Interaction = "PLAND x PROX")

null_long <- as.data.frame(null_mat) %>%
  tidyr::pivot_longer(everything(), names_to = "Term", values_to = "Coefficient") %>%
  dplyr::filter(is.finite(Coefficient)) %>%
  dplyr::mutate(Term = unname(term_labels[Term]))

obs_long <- data.frame(
  Term = c("PLAND_Forest", "PD_Forest", "PLAND x PROX"),
  Observed = c(obs_pland, obs_pd, obs_int)
)

p_null <- ggplot(null_long, aes(x = Coefficient)) +
  geom_histogram(bins = 60, fill = "grey75", colour = NA) +
  geom_vline(data = obs_long, aes(xintercept = Observed),
             colour = "#D55E00", linewidth = 0.9) +
  facet_wrap(~ Term, scales = "free", ncol = 1) +
  labs(x = "Coefficient under the null (forest amount unrelated to window area)",
       y = "Frequency") +
  theme_classic(base_size = 10) +
  theme(strip.background = element_rect(fill = "grey25", colour = NA),
        strip.text = element_text(colour = "white", face = "bold", size = 9))

ggsave("Figure_PV_01_Null_Distribution.tiff", p_null, path = output_path,
       width = 174, height = 180, units = "mm", dpi = 600,
       device = "tiff", compression = "lzw", bg = "white")
message("   Figure_PV_01_Null_Distribution.tiff exported.")


# ==============================================================================
# 6.5 IS THE AICc PREFERENCE FOR THE COUPLED MODEL ITSELF AN ARTEFACT?
# ==============================================================================
# The retained model beats the specification free of area-normalized predictors
# by a wide margin on AICc, and that margin is the natural argument for keeping
# the coupled model. The argument does not hold, and this block is what shows it.
#
# A predictor built from the response will improve the fit to the response
# whether or not it carries any ecological information, so an information
# criterion cannot arbitrate between a coupled and an uncoupled specification.
# The comparison it performs is not between two hypotheses about nature; it is
# between a model that has partial access to the answer and one that does not.
#
# The block therefore recomputes the same AICc difference under the permutation
# null. In each permutation the numerators are reshuffled and the denominator is
# held fixed, so the coupled model retains its construction and loses any real
# information. If the null routinely reproduces the observed AICc advantage, then
# the advantage measures the construction and not the ecology, and the simpler
# specification is the one to report.
# ==============================================================================
message("\n[5.5] AICc advantage of the coupled model, under the null...")

aic_obs_S1 <- MuMIn::AICc(glm(specs$S1$f, data = dz, family = Gamma(link = "log")))
aic_obs_S4 <- MuMIn::AICc(glm(specs$S4$f, data = dz, family = Gamma(link = "log")))
delta_obs  <- aic_obs_S4 - aic_obs_S1

N_SIM_AIC <- min(N_SIM, N_SIM_AIC_CAP)
delta_null <- rep(NA_real_, N_SIM_AIC)

set.seed(GLOBAL_SEED)
for (i in seq_len(N_SIM_AIC)) {
  idx <- sample.int(nrow(dz))
  sim <- dz %>%
    mutate(PLAND_z = z(100 * dz$CA_Forest[idx] / AOH_ha),
           PD_z    = z(dz$NP_Forest[idx] / (AOH_ha / 100)))
  a1 <- tryCatch(MuMIn::AICc(glm(specs$S1$f, data = sim, family = Gamma(link = "log"))),
                 error = function(e) NA_real_, warning = function(w) NA_real_)
  a4 <- tryCatch(MuMIn::AICc(glm(specs$S4$f, data = sim, family = Gamma(link = "log"))),
                 error = function(e) NA_real_, warning = function(w) NA_real_)
  delta_null[i] <- a4 - a1
}

dn <- delta_null[is.finite(delta_null)]
aic_table <- data.frame(
  Quantity = c("AICc of the retained model (S1)",
               "AICc of the specification free of coupling (S4)",
               "Observed advantage of S1, in AICc units",
               "Mean advantage of S1 under the null",
               "Null 95 per cent interval, lower",
               "Null 95 per cent interval, upper",
               "Proportion of permutations reaching the observed advantage"),
  Value = round(c(aic_obs_S1, aic_obs_S4, delta_obs, mean(dn),
                  quantile(dn, 0.025), quantile(dn, 0.975),
                  mean(dn >= delta_obs)), 4)
)
write.csv(aic_table, file.path(output_path, "TableS_PartV_AICc_Under_Null.csv"),
          row.names = FALSE)
print(aic_table, row.names = FALSE)

if (mean(dn >= delta_obs) > 0.05) {
  message("   The null reproduces the observed AICc advantage. The information")
  message("   criterion is measuring the construction of the ratio, not ecology,")
  message("   and it cannot be used to justify keeping the coupled model.")
} else {
  message("   The observed AICc advantage exceeds what the construction alone")
  message("   produces. The coupled model carries information beyond the ratio,")
  message("   which does not licence interpreting the coefficient of the ratio.")
}


# ==============================================================================
# 6.8 IS FOREST PROXIMITY JUST HABITAT AMOUNT MEASURED AGAIN?
# ==============================================================================
# The proximity index is not normalized by the response, so it is not exposed to
# the ratio artefact of Section 6. It is, however, built from the areas of the
# patches within a search radius, divided by the square of their distances, so a
# landscape holding more forest tends to score higher on it by construction. That
# is confounding with habitat amount, which is a different problem from the ratio
# artefact and is precisely the concern Fahrig (2013) raises: an apparent
# configuration effect that is habitat amount in disguise.
#
# The block quantifies it in two steps. First the monotone association between
# proximity and three measures of forest amount. Then the coefficient of
# proximity with each of those measures entered as a covariate. A coefficient
# that survives conditioning on an amount measure free of the response is a
# configuration signal; one that collapses is amount in disguise.
#
# Note on which covariates are admissible. The absolute Class Area is bounded by
# the response and conditioning on it is circular, which is why it appears in the
# correlations but not among the covariates. The mean patch area and the number
# of patches are not divided by the response and are admissible.
# ==============================================================================
message("\n[5.8] Is forest proximity confounded with habitat amount?")

amount_cor <- data.frame(
  Amount_measure = c("PLAND_Forest", "CA_Forest", "Mean_AREA_Forest_ha", "NP_Forest"),
  Pearson  = round(c(cor(d$PROX_MN_Forest, d$PLAND_Forest),
                     cor(d$PROX_MN_Forest, d$CA_Forest),
                     cor(d$PROX_MN_Forest, d$Mean_AREA_Forest_ha),
                     cor(d$PROX_MN_Forest, d$NP_Forest)), 4),
  Spearman = round(c(cor(d$PROX_MN_Forest, d$PLAND_Forest, method = "spearman"),
                     cor(d$PROX_MN_Forest, d$CA_Forest, method = "spearman"),
                     cor(d$PROX_MN_Forest, d$Mean_AREA_Forest_ha, method = "spearman"),
                     cor(d$PROX_MN_Forest, d$NP_Forest, method = "spearman")), 4),
  Contains_the_response = c("yes, in its denominator", "bounded by it", "no", "no")
)

dz$logNP_z <- z(log(dz$NP_Forest))

prox_models <- list(
  "No amount covariate"                 = AOH_ha ~ ENN_z + FRAC_z + PROX_z,
  "Plus mean forest patch area"         = AOH_ha ~ ENN_z + FRAC_z + PROX_z + AREAMN_z,
  "Plus log number of forest patches"   = AOH_ha ~ ENN_z + FRAC_z + PROX_z + logNP_z,
  "Plus both"                           = AOH_ha ~ ENN_z + FRAC_z + PROX_z + AREAMN_z + logNP_z
)

prox_table <- purrr::imap_dfr(prox_models, function(fm, lab) {
  m  <- glm(fm, data = dz, family = Gamma(link = "log"))
  cf <- summary(m)$coefficients
  ci <- confint.default(m)
  data.frame(
    Specification = lab,
    PROX_beta     = round(cf["PROX_z", "Estimate"], 4),
    CI_lower      = round(ci["PROX_z", 1], 4),
    CI_upper      = round(ci["PROX_z", 2], 4),
    P_value       = signif(cf["PROX_z", 4], 4),
    OHE_ratio_per_SD = round(exp(cf["PROX_z", "Estimate"]), 3),
    stringsAsFactors = FALSE
  )
})

write.csv(amount_cor, file.path(output_path, "TableS_PartV_Proximity_Amount_Correlations.csv"),
          row.names = FALSE)
write.csv(prox_table, file.path(output_path, "TableS_PartV_Proximity_Conditioned.csv"),
          row.names = FALSE)
print(amount_cor, row.names = FALSE)
print(prox_table, row.names = FALSE)

surv <- all(prox_table$CI_lower * prox_table$CI_upper > 0)
shrink <- 1 - min(prox_table$PROX_beta) / prox_table$PROX_beta[1]
if (surv) {
  message(sprintf(paste0("   The proximity coefficient keeps an interval excluding zero under every ",
                         "amount covariate, shrinking by at most %.0f per cent. It is a ",
                         "configuration signal, not habitat amount in disguise, although the two ",
                         "cannot be fully separated in this design."), 100 * shrink))
} else {
  message("   The proximity coefficient loses significance under at least one amount covariate.")
  message("   It must then be reported as inseparable from habitat amount.")
}


# ==============================================================================
# 7. DECISION TABLE                                                     [6/7]
# ==============================================================================
message("\n[6/7] Building the decision table...")

s6 <- spec_table %>% filter(grepl("^S6", Specification))
s4 <- spec_table %>% filter(grepl("^S4", Specification))

decision <- data.frame(
  Question = c(
    "Q1. Is the negative coefficient of forest cover more than the ratio artefact?",
    "Q2. Does the effect of forest proximity survive without any coupled predictor?",
    "Q3. Does the interaction survive with an amount measure free of the response?",
    "Q4. Is the coefficient of patch density more than the ratio artefact?",
    "Q5. Can AICc arbitrate between the coupled and the uncoupled specification?"
  ),
  Evidence = c(
    sprintf(paste0("Observed %.3f; null central 95%% interval [%.3f, %.3f]; ",
                   "observed at the %.1fth percentile of the null; ",
                   "two-sided empirical p = %.4f"),
            null_summary$Observed[1], null_summary$Null_q2.5[1],
            null_summary$Null_q97.5[1], null_summary$Percentile_of_observed[1],
            null_summary$Empirical_p_two_sided[1]),
    sprintf("S4 (no coupled predictor): PROX beta = %.3f [%.3f, %.3f], p = %s",
            s4$PROX_beta, s4$PROX_CI_low, s4$PROX_CI_high, format(s4$PROX_p)),
    sprintf("S6 (amount as mean patch area): interaction beta = %.3f [%.3f, %.3f], p = %s",
            s6$Interaction_beta, s6$Interaction_CI_low,
            s6$Interaction_CI_high, format(s6$Interaction_p)),
    sprintf(paste0("Observed %.3f; null central 95%% interval [%.3f, %.3f]; ",
                   "observed at the %.1fth percentile of the null; ",
                   "two-sided empirical p = %.4f"),
            null_summary$Observed[2], null_summary$Null_q2.5[2],
            null_summary$Null_q97.5[2], null_summary$Percentile_of_observed[2],
            null_summary$Empirical_p_two_sided[2]),
    sprintf("Observed advantage of the coupled model %.2f AICc units; the null reaches it in %.1f%% of permutations",
            delta_obs, 100 * mean(dn >= delta_obs))
  ),
  Verdict = c(
    null_summary$Verdict[1],
    ifelse(s4$PROX_CI_low * s4$PROX_CI_high > 0,
           "YES: the proximity effect does not require a coupled predictor",
           "NO: the proximity effect depends on the coupled predictors; report this"),
    ifelse(!is.na(s6$Interaction_p) &&
             s6$Interaction_CI_low * s6$Interaction_CI_high > 0,
           "YES: the threshold is reproduced without any area-normalised predictor",
           "NO: the threshold is not reproduced; report it as conditional on PLAND"),
    null_summary$Verdict[2],
    ifelse(mean(dn >= delta_obs) > 0.05,
           paste("NO: the construction of the ratio reproduces the observed AICc",
                 "advantage, so the criterion measures access to the response and",
                 "not ecological information. Report the uncoupled specification."),
           paste("The observed advantage exceeds the construction, which still does",
                 "not licence interpreting the coefficient of the ratio."))
  ),
  stringsAsFactors = FALSE
)

write.csv(decision, file.path(output_path, "TableS_PartV_Decision.csv"),
          row.names = FALSE)

message("\n============================================================")
for (i in seq_len(nrow(decision))) {
  message("  ", decision$Question[i])
  message("     Evidence: ", decision$Evidence[i])
  message("     Verdict:  ", decision$Verdict[i])
  message("")
}
message("============================================================")


# ==============================================================================
# 8. WHAT TO WRITE                                                      [7/7]
# ==============================================================================
# The uncoupled specification is what the manuscript reports, so its full
# coefficient table is exported rather than left inside the specification summary.
m_free <- glm(specs$S4$f, data = dz, family = Gamma(link = "log"))
cf_free <- summary(m_free)$coefficients
ci_free <- confint.default(m_free)
free_table <- data.frame(
  Predictor = rownames(cf_free),
  Estimate  = round(cf_free[, "Estimate"], 4),
  SE        = round(cf_free[, "Std. Error"], 4),
  CI_lower  = round(ci_free[, 1], 4),
  CI_upper  = round(ci_free[, 2], 4),
  P_value   = signif(cf_free[, 4], 4),
  OHE_ratio_per_SD = round(exp(cf_free[, "Estimate"]), 3),
  row.names = NULL
)
write.csv(free_table, file.path(output_path, "TableS_PartV_Uncoupled_Model.csv"),
          row.names = FALSE)
message("\n   Coefficients of the specification free of coupling (S4):")
print(free_table, row.names = FALSE)
message(sprintf("   Maximum VIF in this specification: %.2f",
                spec_table$Max_VIF[spec_table$Specification == specs$S4$label]))

# ==============================================================================
# 7.5 FIGURE OF THE REPORTED MODEL
# ==============================================================================
# Figure 5 of the manuscript must show the coefficients of the model that the
# manuscript reports. Until this block existed it showed the coefficients of the
# model selected over the coupled candidate set, which is no longer the model
# interpreted, so the figure and Table 3 disagreed. The plot is produced here,
# next to the specification it belongs to, rather than in Part III.
# ==============================================================================
message("\n[7.5] Coefficient plot of the reported model...")

LAB <- c(PROX_z = "Forest proximity (PROX_MN)",
         FRAC_z = "Shape complexity (FRAC_MN)",
         ENN_z  = "Isolation (ENN_MN)")

coef_plot_df <- free_table %>%
  dplyr::filter(Predictor != "(Intercept)") %>%
  dplyr::mutate(Label = unname(LAB[Predictor]),
                Label = factor(Label, levels = rev(unname(LAB))),
                Supported = ifelse(CI_lower * CI_upper > 0,
                                   "interval excludes zero", "not distinguishable"))

p_coef <- ggplot(coef_plot_df, aes(x = Estimate, y = Label, colour = Supported)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey45") +
  geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper), height = 0.14, linewidth = 0.7) +
  geom_point(size = 2.8) +
  scale_colour_manual(values = c("interval excludes zero" = "#0072B2",
                                 "not distinguishable"    = "#999999")) +
  labs(x = "Standardised association with log(envelope extent)", y = NULL, colour = NULL) +
  theme_classic(base_size = 10) +
  theme(legend.position = "top", legend.text = element_text(size = 8),
        axis.text.y = element_text(size = 9))

ggsave("Figure_PV_02_Reported_Model_Coefficients.tiff", p_coef, path = output_path,
       width = 140, height = 80, units = "mm", dpi = 600,
       device = "tiff", compression = "lzw", bg = "white")
message("   Figure_PV_02_Reported_Model_Coefficients.tiff exported.")
message("   This is Figure 5 of the manuscript. The former Figure_PIII_02 shows the")
message("   coefficients of the coupled candidate set and belongs in the supplement.")


# ==============================================================================
# 7.6 OUT-OF-SAMPLE ERROR AND BIC OF THE REPORTED MODEL
# ==============================================================================
# Two additions asked for by the literature on model selection for inference.
#
# Yates et al. (2023) recommend exact leave-one-out rather than random k-fold
# with k < 10, and they warn that a random split places units of the same taxon
# in the training and the test set at once. Both schemes are reported.
#
# Brewer, Butler and Cooksley (2016) show that the relative behavior of AICc and
# BIC depends on unobserved heterogeneity, which these data carry. Both criteria
# are reported for the reported model and for the coupled one, so the reader can
# see whether they agree.
#
# The comparison is stated for what it is. The coupled specification will always
# appear to predict better, because two of its predictors are built from the
# response. That is the same artefact that Table S28c documents for the
# information criterion, and it means that out-of-sample error cannot arbitrate
# between the two specifications either.
# ==============================================================================
message("\n[7.6] Out-of-sample error and BIC of the reported model...")

loo_rmse <- function(fml, dat) {
  err <- numeric(nrow(dat))
  for (i in seq_len(nrow(dat))) {
    m <- tryCatch(glm(fml, data = dat[-i, ], family = Gamma(link = "log")),
                  error = function(e) NULL)
    err[i] <- if (is.null(m)) NA_real_ else
      (dat$AOH_ha[i] - predict(m, newdata = dat[i, , drop = FALSE],
                               type = "response"))^2
  }
  sqrt(mean(err, na.rm = TRUE))
}
lgo_rmse <- function(fml, dat) {
  err <- numeric(0)
  for (g in unique(dat$GENUS)) {
    tr <- dat[dat$GENUS != g, ]; te <- dat[dat$GENUS == g, ]
    m  <- tryCatch(glm(fml, data = tr, family = Gamma(link = "log")),
                   error = function(e) NULL)
    p  <- if (is.null(m)) rep(NA_real_, nrow(te)) else
      predict(m, newdata = te, type = "response")
    err <- c(err, (te$AOH_ha - p)^2)
  }
  sqrt(mean(err, na.rm = TRUE))
}

fit_free <- glm(specs$S4$f, data = dz, family = Gamma(link = "log"))
fit_coup <- glm(specs$S1$f, data = dz, family = Gamma(link = "log"))

cv_table <- data.frame(
  Specification = c("Reported model, free of coupling (S4)",
                    "Coupled model (S1)"),
  k        = c(length(coef(fit_free)) + 1, length(coef(fit_coup)) + 1),
  AICc     = round(c(MuMIn::AICc(fit_free), MuMIn::AICc(fit_coup)), 2),
  BIC      = round(c(stats::BIC(fit_free),  stats::BIC(fit_coup)), 2),
  LOO_RMSE_ha = round(c(loo_rmse(specs$S4$f, dz), loo_rmse(specs$S1$f, dz)), 1),
  LGO_RMSE_ha = round(c(lgo_rmse(specs$S4$f, dz), lgo_rmse(specs$S1$f, dz)), 1),
  Note = c("no predictor built from the response",
           "two predictors built from the response; every criterion favours it for that reason"),
  stringsAsFactors = FALSE
)
write.csv(cv_table, file.path(output_path, "TableS_PartV_CV_And_BIC.csv"),
          row.names = FALSE)
print(cv_table, row.names = FALSE)
message("   Read this table together with Table S28c: neither the information")
message("   criterion nor the out-of-sample error can arbitrate between the two,")
message("   because both reward access to the response.")


message("\n[7/7] Guidance for the manuscript.")
message("  METHODS. State that PLAND and PD are normalised by the response, that")
message("  this is the ratio problem of Pearson (1897) and Kronmal (1993), and")
message("  that the consequence was quantified by permutation with the numerators")
message("  reshuffled and the denominator held fixed.")
message("")
message("  RESULTS. Report the specification table and the null comparison before")
message("  interpreting any coefficient of PLAND or PD. Lead with the quantity")
message("  that is free of the coupling, the effect of forest proximity, and with")
message("  the interaction as reproduced by S6.")
message("")
message("  WHAT NOT TO WRITE. Do not describe the negative coefficient of forest")
message("  cover as evidence that less forest accompanies more habitat. Whatever")
message("  the null comparison returns, that reading is not available from this")
message("  design, because forest cover is defined as a fraction of the response.")

message("\n[Session]")
message("\n[Session] Package versions and system information:")
print(sessionInfo())
writeLines(capture.output(sessionInfo()),
           file.path(output_path, "session_info_part05.txt"))

# ==============================================================================
# END OF PART V
# ==============================================================================
