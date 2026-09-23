# What each script does, and why

Manuscript: *Assessing the Landscape Composition and Configuration of Threatened Medium- and Large-Sized Mammals using Habitat Envelopes as Sampling Windows*

One section per script. Each states the question the script answers, the decisions it encodes, the
guards it enforces, and what it writes.

---

## `00_recompute_class_metrics.R`

**Question.** Do the class-level metrics stored in the consolidated dataset reproduce from the primary
rasters, and can the empty cells of two sampling units be filled without changing method?

**Why it exists.** Two units, *A. guariba* UA 6 and *B. arachnoides* UA 1, reached the consolidated
dataset with empty cells in edge density and largest patch index of every class, and in number of
patches, patch density and proximity index of the farming class. They had been inserted by a
reconstruction that could only read the exported FRAGSTATS tables, and those tables carry neither patch
perimeter nor patch counts for the union of pasture and agriculture. Downstream, `complete.cases()`
dropped both units without a message, which is where the sample sizes of 63 and 65 in earlier versions
came from.

**What it does.** Recomputes every class-level metric for all 67 units with `landscapemetrics`, using
the same reclassification and the same eight-neighbour rule that produced the original values. Recovers
the proximity index from the FRAGSTATS patch tables by the same pooled mean over pasture and
agriculture patches. Then compares, metric by metric, against the stored values.

**Nothing is imputed.** Every filled value is measured again from the same raster.

**One performance note that matters.** The per-unit rasters cover the bounding box of every unit of the
taxon, so more than 99 per cent of each file is NoData padding: 5.33 billion pixels across the 67 files,
of which about 0.7 per cent carry data. The script crops each raster with `terra::trim()` before any
metric is computed, which reduces the volume by a factor of 20 to 915, and adds a one-cell NoData margin
back with `terra::extend()` so that edge density cannot react to the crop. The equality is then verified
against the stored values rather than assumed.

**Guards.** Stops unless 67 rasters are found. Stops on a duplicated (taxon, unit) key. Caches each
unit, so an interrupted run resumes.

**Writes.** `Tabela_Metricas_Recomputada.csv`, `TableS_Recompute_Agreement.csv`,
`TableS_Repaired_Cells.csv`.

---

## `01_prepare_data_and_ordination.R`

**Question.** What is the composition of the landscape inside each occupied habitat envelope, and does
it differ among taxonomic and functional groups?

**What it does, in order.**

1. **Canonicalises taxon labels.** `canonical_taxon()` maps every orthographic variant found in the
   project onto one label. The variant that mattered was *Brachyteles hypoxanthus*: the tabular data
   carried `B_hypoxantus` and the shapefile carried the valid spelling, so the compound key never
   matched and the two units of that species reached the spatial analyses without coordinates. An
   unrecognised name raises an error rather than passing through.
2. **Repairs the metric gaps** from the table written by script 00.
3. **Recodes structural zeros, conditionally.** FRAGSTATS omits the row of a class absent from a unit,
   and the consolidation turned that omission into `NA`. A unit with no herbaceous patch has zero class
   area, not an unknown value. But a unit with 91.6 per cent forest cover and an empty edge density cell
   has a missing value, not a zero. The previous pipeline recoded every `NA` to zero and thereby wrote
   two false zeros of forest edge density into the variation partitioning. A cell is now zeroed only
   when the class area of its class is zero.
4. **Extracts polygon centroids** from the per-species shapefiles, repairing invalid geometry and
   falling back to `st_point_on_surface()` where a centroid falls outside its polygon.
5. **Ordinates and tests.** Non-metric multidimensional scaling on Bray-Curtis dissimilarities of
   square-root transformed class area; `envfit` for the configuration metrics; `betadisper` with a
   permutation test before the permutational analysis of variance, because the design is unbalanced and
   the test responds to dispersion as well as to location.

**Guards.** Stops unless the input has 67 rows, unless all 67 carry coordinates, and unless all 67 are
complete over the fourteen-predictor analysis pool used by scripts 02, 03 and 05.

**Writes.** `Data_Raw_WithCoords.csv`, which is the sole input of scripts 02, 03 and 05. Tables 1 to
4, Tables S2 to S5, `TableS_NMDS_Summary.csv`, `TableS_Metric_Repair_Log.csv`,
`TableS_Sample_Sizes_PartI.csv`, and Figures 1 to 3.

---

## `02_constrained_ordination_and_partitioning.R`

**Question.** Which landscape and environmental predictors structure the composition, and how much of
the variation in envelope extent is attributable to the landscape once spatial structure is accounted
for?

**What it does.**

- **Redundancy analysis** of the square-root transformed class area of the three focal classes,
  constrained by configuration and environmental predictors. The percentage of landscape is excluded
  from this pool because it is a transformation of the ordination response.
- **Variation partitioning** of log(envelope extent) between a landscape set and a spatial set. The
  spatial set is built from distance-based Moran eigenvector maps of the unit centroids, modelling
  positive autocorrelation only. Individual eigenvectors are selected only if the global spatial model
  is significant, which controls the Type I error inflation that arises otherwise.
- **Iterative variance inflation reduction** at a threshold of 5, following Zuur et al. (2010). Two
  alternative strategies, principal components and forward selection, are present as documented but
  inactive blocks.

**What it no longer does.** The emergency centroid rebuild was removed. It duplicated the logic of
script 01 in a second place and the two copies drifted apart.

**Guards.** Verifies the input rather than repairing it. Stops if any unit is incomplete over the
analysis pool, if any lacks coordinates, or if edge density is zero in a unit where the class is
present, which is the signature of the unconditional structural zero.

---

## `03_univariate_model_selection.R`

**Question.** Which forest configuration metrics accompany a larger envelope, under the candidate set
defined a priori?

**What it does.**

- Compares six specifications on a common response scale: ordinary least squares, a linear mixed model
  with a random intercept for the genus, and a Gamma generalized linear model with a logarithmic link,
  each with and without the a priori interaction. Compared by AICc, by BIC, and by exact leave-one-out
  cross validation with leave-one-genus-out reported beside it. Both criteria appear because their
  relative behaviour depends on unobserved heterogeneity, which 67 units from 9 taxa certainly carry
  (Brewer et al., 2016). Leave-one-out replaces the earlier random ten-fold scheme, which had neither
  the low bias of leave-one-out nor a bias correction (Yates et al., 2023).
- **Fits the full a priori model once, without selection, and reports its coefficients as the
  inferential result.** Yates et al. (2023) show that inference on coefficients is biased when it
  follows a selection step, and that valid inference needs a single carefully specified model.
  Tredennick et al. (2021) make the same point by goal. The candidate set here was pre-specified, so
  the full model is that single specification.
- **Refits the confirmatory model with the Moran eigenvectors retained by script 02**, and compares the
  coefficients. Diniz-Filho et al. (2008) show that an information criterion computed on a non-spatial
  regression is sensitive to spatial autocorrelation and yields unstable minimum adequate models. The
  selection here is non-spatial and the data carry detectable spatial structure, so the check is run
  rather than left as a caveat.
- Explores the candidate set exhaustively with `dredge`, and reports both the conditional and the full
  average, with the revised variance estimator, alongside the sum of Akaike weights. This block is
  labelled as a description of selection uncertainty, not as a second set of estimates to test.
- Screens the dependence among units of the same taxon by comparing against the mixed specification and
  reporting the intraclass correlation alongside the Akaike weight.
- **Estimates the genus-level slopes twice, by two designs that disagree on purpose.** Section 10.6 fits
  a univariate regression inside each genus, so a genus of four units estimates its own residual
  variance from four units. Section 10.6b fits one Gamma model per metric to all 67 units with the genus
  as a moderator, reads each genus slope as the linear contrast of the reference slope and its own
  interaction term, and takes the standard error from the full covariance matrix rather than the
  diagonal. The second design lets the small genera borrow precision from the rest, and the two bracket
  the result. Where they disagree, the disagreement is about how much a small genus may speak for
  itself, and both are reported for that reason. Section 10.6b uses only base R; it replaces the
  marginal-effects script that was removed, and without it Table S29 of the supplement had no producing
  code.
- Checks residual behaviour by simulation with `DHARMa`.

**The a priori interaction.** The fragmentation threshold hypothesis is a statement about an interaction
between connectivity and habitat amount, not about either term alone, so the product of the proximity
index and the percentage of forest is admitted. It is the only interaction considered.

**Three outputs, three roles.** `Table_S_PartIII_ConfirmatoryModel.csv` carries the inferential
estimates. `Table_S_PartIII_SpatialSensitivity.csv` shows whether they hold once spatial structure is
in the model. `Table_S_PartIII_Averaging_FullVsConditional.csv` describes how much the selection was
uncertain. Only the first is inference.

**Two more, for the genus-level result.** `Table_S_PartIII_EffectByGenus.csv` holds the within-genus
slopes and `Table_S_PartIII_SlopeByGenus_Moderator.csv` holds the moderated ones. Per-genus samples
range from four to twelve units, so neither is confirmatory.

**Figure labels fail loudly.** `predictor_labels` covers every predictor, including forest cover and
the interaction term, and `label_predictors()` returns the variable name with a warning when a label is
missing rather than returning `NA`. `strip_z()` removes the standardisation suffix from every term of
an interaction, not only from the last one. A `stopifnot` checks the heatmap labels before any figure
is drawn. Together these replace a silent failure that printed an axis reading `NA`.

**Important.** The model selected here is **not** the model the manuscript interprets. Script 05
explains why. This script is retained because the candidate set was pre-specified and reporting the
pre-specified analysis is part of the record.

---

## `04_geometric_coupling_diagnostic.R`

**Question.** Which predictors are algebraically tied to the response, and by how much?

**What it does.** Verifies numerically, on all 67 units, the identity that links patch density to the
response: log(extent) equals log(number of patches) plus log(100) minus log(patch density). Then fits
four models: the reference, a legitimate control on the percentage of landscape, and two circular
controls on the number of patches and on the class area, the last two reported precisely to make the
circularity visible.

---

## `05_ratio_artefact_test.R`

**Question.** Is the negative coefficient of forest cover more than the ratio structure produces on its
own, and does the fragmentation threshold survive without the coupled predictors?

**This script decides which model the manuscript reports.**

**The permutation null.** The class area and the number of forest patches are reshuffled jointly across
sampling units, preserving their mutual association, while the area of the unit is held at its observed
value. The percentage of landscape and the patch density are recomputed from the permuted numerators and
the observed denominator, and the model is refitted, 9,999 times with a fixed seed. Under this null
there is by construction no relationship between the amount of forest and the extent of the envelope,
yet the ratio structure remains, so any coefficient obtained measures the construction alone.

**The same null is applied to the information criterion.** This is the part most easily overlooked.
Model selection is routinely treated as a neutral arbiter, and here it is not: a predictor built from the
response improves the fit to the response whether or not it carries ecological information.

**Six alternative specifications** remove the normalised metrics singly and jointly and measure habitat
amount instead by the absolute class area and by the mean forest patch area.

**The interaction is always proximity times habitat amount.** What changes across the six is which
variable measures amount: the percentage of landscape in S1 and S3, which carries the response in its
denominator; the log class area in S5, which the response bounds from above; and the mean forest patch
area in S6, the only one the response neither contains nor bounds. S2 and S4 carry no amount term, so
no product can be formed in them. Shape complexity never enters an interaction; it is a main effect in
all six.

**One further comparison, added because script 03 cannot make it.** Script 03 finds the interaction
structure better than the additive one by more than thirteen AICc units, but it makes that comparison
with the percentage of landscape in both arms. Section 3b repeats the comparison twice, once with the
coupled measure of amount and once with the mean patch area, and writes
`TableS_PartV_Interaction_Gain.csv`. The gap between the two differences is the part of the advantage
that the normalisation produces rather than any ecological threshold.

**A separate screen for the proximity index.** The index is not divided by the area of the window, so it
is not exposed to the ratio artefact, but it is built from the areas of neighbouring patches, so a
landscape holding more forest scores higher on it by construction. The script reports the association
with four measures of forest amount and refits the model with those that are free of the response.

**Writes,** among others, the coefficient table and the coefficient plot of the model the manuscript
reports.

---

## `08_sampling_effort_sensitivity.R`

The response is the area of a kernel utilisation distribution at the 99% isopleth with a single
smoothing parameter. At fixed h that area grows with the number and the spread of the records that
produced it, and effort in secondary records concentrates near roads, research stations and protected
areas, which are also where forest persists. Part of the association between configuration and envelope
extent could therefore be effort rather than landscape structure. Appendix S1 addresses the spatial
unevenness of effort; it does not address this dependence, and this script measures it.

The measure is the number of retained records falling inside each unit, entered as `log(records + 1)`
and standardised like the other predictors. The model of Table 3 is refitted with that covariate added.

**How to read the result.** The covariate is not a pure nuisance. Section 2.3 of the manuscript takes a
larger, more stable population to maintain a broader area of use, and such a population also generates
more records, so the count partly mediates the same latent quantity as the response. Conditioning on it
removes effort and part of the ecological signal at once. The augmented coefficients are therefore a
lower bound and those of Table 3 an upper bound. Both associations keep their sign and exclude zero at
both bounds.

**What the measure is not.** Only the Rio de Janeiro subset of the occurrence compilation is archived
here, while the kernels were estimated on the biome-wide compilation of Macedo et al. (2019). The count
is the density of records inside each envelope, not the complete input to the kernel. Ten units hold no
record of the state subset, which is why `log(records + 1)` is used rather than `log(records)`, and why
the analysis is repeated on the 57 units that hold at least one.

**Writes** `PartVI/TableS42_Sampling_Effort.csv`, `PartVI/TableS42b_Sampling_Effort_Robustness.csv` and
`PartVI/Records_Per_Unit.csv`. It reads the model of Table 3 and therefore runs after script 05.

---

## `99_legacy_metric_extraction.R`

The original extraction, superseded by script 00. Kept because it documents the provenance of the
consolidated dataset and the class aggregation, and because the proximity index is still recovered from
FRAGSTATS by the procedure it describes. It should not be run as part of the current chain.

---

## Why there is no marginal-effects script

An earlier version of the chain contained a script that computed simple slopes with `emmeans` and
adjusted predictions with `ggeffects`, together with a Johnson-Neyman interval. It was removed.

The reason is that the model the manuscript reports contains no interaction, so its coefficients are
already the associations, and a marginal effect adds nothing to them. The script existed to interpret
the interaction of the a priori candidate set, and script 05 shows that the interaction is not
separable from the geometry of the design. Reporting the marginal effects of a term that is not
interpreted would have given it a prominence the analysis does not support.

The simple slopes and the Johnson-Neyman boundary that the script produced are therefore no longer
reported, in the manuscript or in the supplement.

---

## What has no script here

The weighted conservation value index, which supports Figure 7 of the manuscript and Tables S36 and
S37 of the supplement, was built in ArcGIS Pro 3.1.0 rather than in this chain. The steps are stated in
Section 5 of the supplement: the nine envelope polygons are intersected under the South America Albers
Equal Area Conic projection, the sensitivity weights are summed over the resulting polygons, the index
surface is classified, and the result is intersected with the protected area network of the Instituto
Estadual do Meio Ambiente. The input polygons are the ones this chain uses and are archived with the
repository, so the index can be rebuilt from them.

Figure 1, the study area map, was also produced in a geographic information system and has no script.

---

## Guards common to the chain

Every script stops rather than analysing a reduced sample. The specific checks are:

- the input has 67 rows;
- all 67 carry centroid coordinates;
- all 67 are complete over the analysis pool;
- no metric was zeroed in a class that is present in the unit;
- the sample size of every analysis equals 67, verified and written to disk.

The reason for this severity is historical. Three separate defects each removed units silently, and each
surfaced only later as a reduced sample size in a supplementary table.

---

## The connectivity demonstration is not archived here

An exploratory calculation of the probability of connectivity of Saura and Pascual-Hortal (2007), and of
the importance of each forest patch within an envelope, dPC, was prepared on 25 August 2026 as material
for the defence. It is reported in neither the manuscript nor the supplement, and **the scripts and their
outputs are not archived in this repository**. This section records the argument, because it bears
directly on Section 2.5.1 of the manuscript, and not the code.

PC is normalised by the total landscape area, which in this design is the area of the envelope and
therefore the response, so PC must not be compared across sampling units. dPC is the relative drop in PC
when one patch is removed, and the squared landscape area cancels between numerator and denominator, so
removing a patch does not change it. The ranking of patches inside one envelope is therefore free of the
geometric coupling the manuscript describes, while a ranking of whole units by PC would not be.

Anyone rebuilding that demonstration should note that the dispersal distance is assumed rather than
measured, so the ranking has to be reported at more than one median distance, and that the connector
fraction of Saura and Rubio (2010) is identically zero under a complete graph with straight-line
distances and a negative exponential kernel, which is why links have to be cut above a threshold.
