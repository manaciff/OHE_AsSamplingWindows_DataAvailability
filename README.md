# Assessing the Landscape Composition and Configuration of Threatened Medium- and Large-Sized Mammals using Habitat Envelopes as Sampling Windows

Short title: Habitat envelopes as sampling windows.

Code and data for the analysis of landscape configuration and the extent of the occupied habitat
envelope of nine threatened medium and large mammal taxa in the Atlantic Forest of Rio de Janeiro
State, Brazil.

**Author:** Maria Eduarda Nacif
**Status:** manuscript in preparation. Last revision of the analysis chain: 22 August 2026.
**Language:** R 4.5.2

---

## What the study asks

Two questions, in this order.

**The goal is inference, not prediction.** That is stated here because it governs the whole workflow.
Tredennick et al. (2021) show that exploration, inference and prediction require different model
selection procedures, and that most confusion comes from not declaring the goal first. The estimates
this repository reports come from single models specified a priori and fitted once. The exhaustive
search over subsets is present, and is labelled as a description of selection uncertainty rather than as
a second set of estimates to test, because inference on coefficients is biased when it follows selection
(Yates et al., 2023). These data cannot support prediction: every predictor is measured inside the
polygon whose area is the response, so no predictor exists before the response is known.

**Methodological.** Landscape studies often measure habitat configuration inside a window whose area is
also the response variable. Home ranges, utilisation distributions, territories and occupied polygons
are all used this way. When that happens, any metric normalised by the area of the window carries the
response in its own denominator. This repository quantifies the consequence and shows that an
information criterion cannot arbitrate between a model that contains such a metric and one that does
not.

**Empirical.** Which features of the remaining forest accompany a larger occupied habitat envelope, once
the metrics exposed to that dependence are set aside?

## The response variable, and what it is not

The response is the **area of the occupied habitat envelope (OHE)**: the 99 per cent isopleth of a
kernel utilisation distribution fitted to occurrence records, masked and split into disjoint parts, one
per local population. It is measured in hectares.

It is **not** the Area of Habitat of Brooks et al. (2019), which is the suitable habitat within a
species range. In these data the forest inside the envelope is the closer analogue of that quantity and
occupies a median of 54.9 per cent of the envelope. It is **not** the Area of Occupancy of the IUCN Red
List, which is measured on a mandatory 2 by 2 km grid. No IUCN threshold is applied to any value here.

## Sampling design

67 sampling units, from 9 taxa in 8 genera. Every analysis in this repository uses all 67; each script
stops with a named list of offending cells if any step would reduce that number.

---

## Repository layout

The scripts resolve every path with `here()` from the root of the repository, so the folder names
below are the ones the code expects. A `.here` file marks the root for anyone who downloads a ZIP
rather than cloning.

```
R/                                 analysis scripts, numbered in execution order
run_all.R                          runs the chain in order; scripts 00 and 01 are commented out
                                   and run separately, see the note inside the file
Dados/Processados/                 tabular inputs, and the files the scripts write back
Dados/FRAGSTATS_RESULT/            FRAGSTATS patch and class tables, read by script 00
Dados/Rasters/UA_RASTER/           67 per-unit rasters, read by script 00; see the note inside
Dados/Shapefiles/                  hexagonal grid and the conservation value index layers
08_Dados_Especies/                 per-taxon envelope polygons, read by script 01 for the centroids
Outputs/Manuscrito/                tables and figures, one subfolder per part of the chain
docs/                              per-script documentation and the table map
arcpy_AOH_area_full_precision.py   exports the response at full precision from ArcGIS Pro
```

Everything except the rasters is tracked here, about 29 MB in total. Scripts 01 to 05 run on the
tabular files alone, and they produce every table and figure of the manuscript and of the
supplementary material. Script 00 additionally needs the 2.6 GB of rasters, which are archived
separately; `Dados/Rasters/UA_RASTER/READ_THIS_FIRST.md` says where to get them and why they are
not here.

## Execution order

The scripts must run in numerical order. Each one verifies what it receives from the previous one and
stops rather than analysing a reduced sample.

| # | Script | Reads | Writes | Runtime |
|---|---|---|---|---|
| 00 | `00_recompute_class_metrics.R` | 67 rasters, FRAGSTATS patch tables | `Tabela_Metricas_Recomputada.csv`, agreement table | minutes |
| 01 | `01_prepare_data_and_ordination.R` | raw table, recomputed table, shapefiles | `Data_Raw_WithCoords.csv`, Tables 1 to 4, Figures 1 to 3 | seconds |
| 02 | `02_constrained_ordination_and_partitioning.R` | `Data_Raw_WithCoords.csv` | RDA, variation partitioning, Moran eigenvectors | ~1 min |
| 03 | `03_univariate_model_selection.R` | `Data_Raw_WithCoords.csv`, MEM vectors from 02 | confirmatory model, spatial sensitivity, genus-level slopes, selection uncertainty | ~3 min |
| 04 | `04_geometric_coupling_diagnostic.R` | `Data_Raw_FINAL.csv` | algebraic identity, circular controls | seconds |
| 05 | `05_ratio_artefact_test.R` | `Data_Raw_WithCoords.csv` | permutation null, reported model and its figure | ~10 min |
| 99 | `99_legacy_metric_extraction.R` | rasters | superseded by 00; kept for provenance | — |

**Run 02 before 03.** Script 02 writes the Moran eigenvectors that the spatial sensitivity check of
script 03 reads.

**Do not skip 00.** Two sampling units reach the consolidated dataset with empty cells, and without 00
they are dropped silently by `complete.cases()` further downstream.

## Reproducibility

- Seed fixed at 123 in every script.
- 9,999 permutations in every permutation test.
- Cross validation is exact leave-one-out, with leave-one-genus-out reported beside it.
- AICc and BIC are both reported; where they disagree the disagreement is the result.
- Paths are relative, resolved with `here::here()`.
- `sessionInfo()` is printed at the end of every script and written to `outputs/session_info_part*.txt`.
- Term tests in the constrained ordination are marginal, not sequential.
- The ordination response is Hellinger transformed, so that the size of the sampling window does not
  enter an analysis of composition.

### Packages

`here`, `terra`, `landscapemetrics`, `sf`, `tidyverse`, `vegan`, `adespatial`, `car`, `MuMIn`,
`DHARMa`, `lme4`, `ggridges`, `ggtext`, `ggrepel`, `corrplot`, `patchwork`,
`scales`.

### Two packages that must not be attached carelessly

`MASS` exports `select()` and `car` exports `recode()`, both of which mask the `dplyr` functions of the
same name. Script 05 avoids `ppcor` for that reason and qualifies every call as `dplyr::select()`.

## Data sources

- **Land cover:** MapBiomas Collection 10.1, reference year 2024, 30 m, minimum mapping unit near 0.5
  ha. Reported accuracy for the Atlantic Forest is 91.5 per cent at level 1 and 86.1 per cent at level 2.
- **Occurrence records:** compiled by Macedo et al. (2019), filtered to 1990 to 2020.
- **Landscape metrics:** `landscapemetrics` for class-level metrics; FRAGSTATS for the proximity index,
  which requires a user-defined search radius that `landscapemetrics` does not implement.

## Known limitations of the design

Stated here because they govern what the code can and cannot demonstrate.

1. The sampling window is the response. Metrics normalised by it are not separable from that
   construction. Script 05 measures how far this reaches.
2. The occurrence records span 1990 to 2020 while the land cover is 2024.
3. A single kernel bandwidth, calibrated on *Panthera onca*, is applied to all taxa.
4. Per-genus samples range from four to twelve units, so per-genus estimates are exploratory.

## What is not in this repository

The weighted conservation value index behind Figure 8 and Tables S36 and S37, and the study area map of
Figure 1, were built in a geographic information system rather than in R. Section 5 of the supplement
states the geoprocessing steps, and the input polygons are archived here, so both can be rebuilt.

## What must be in place before this repository is archived

Two things are still outstanding.

1. **The 67 per-unit rasters**, 2.6 GB, which only script 00 reads. Deposit them on Zenodo and record
   the DOI in `CITATION.cff`. Scripts 01 to 05 do not need them, and the result of the recomputation
   script 00 performs is already tracked in `Outputs/Manuscrito/Part0/`.
2. **A completed `CITATION.cff`**: the release date, the repository URL and the archive DOI are still
   placeholders, and the same three values belong in the Data Availability statement of the
   manuscript.

## Revision of 24 August 2026

- **The envelope protocol is attributed.** Section 2.4 of Macedo et al. (2019) describes the same
  procedure this study executed, down to the reference species and the bandwidth rule. The manuscript
  cited that paper only as the source of the occurrence records; Section 2.3 now states the filiation
  and marks the one departure, the 99% isopleth in place of the whole utilization distribution.
- **Figure S10.** The coefficient plot that script `03` had always written was reproduced in neither
  document. It is now Figure S10, redrawn with the three terms the response normalises separated from
  the three that can be read, and without a significance mark on the first three.
- **Eight defects in the chain**, from a guard that fired with the wrong error to a manifest that
  advertised a file no line writes. `docs/Registro_de_Correcoes_24.08.md` lists each one with its
  effect.
- **Dead code removed** and the six scripts standardised: one export of session information per
  script, constants in the constants block, no unreachable branches, no packages attached without
  use. Script `99` can no longer be run by accident over the inputs of the chain.
- Script `04`, previously bilingual, is in English throughout. It was re-run on the same input and
  reproduces the verified output value for value.

## Revision of 23 August 2026

The response is now read at full precision everywhere, and the chain was re-run from script 00. The
changes it produced are small and no conclusion moved; `docs/Registro_de_Correcoes_23.08.md` lists them
one by one.

- `arcpy_AOH_area_full_precision.py` exports the area of every sampling unit from the same polygons,
  with every significant digit. Script 01 consumes the file when it is present, joins it by
  `(SPECIES, UA_ID)`, and stops if any unit is missing or departs by more than five per cent.
- **Script 04 read a different file from scripts 03 and 05.** It read `Data_Raw_FINAL.csv`, in which the
  response is stored as a whole number, so the two tables that verify the same algebraic identity
  disagreed: Table S26 reported a maximum discrepancy of 0.0290 and Table S28e reported 0.0277. Script
  04 now reads `Data_Raw_WithCoords.csv`, the file the rest of the chain reads, and both tables report
  0.0277.
- **The captions of Figures S5 and S6 described a mixed model.** Section 8 of script 03 retains the
  Gamma model on these data, so `global_model` is a `glm` and both figures are panels of it. Script 03
  now exports `Table_S_PartIII_ObsPred_R2.csv`, so the values quoted in the caption of Figure S5 come
  from a file rather than from the image.
- **`run_all.R` could not find the scripts.** With `ANALISES_TCC.Rproj` open, `here()` resolves to the
  working folder, where `R/` does not exist, and `source(here("R", s))` failed on the first script.
  `run_all.R` now locates the chain from its own position and checks that `Dados/` and `Outputs/` exist
  before the first `source`.
- The exported row labels and the figure axis labels now call the response OHE, the name the
  manuscript uses, instead of AOH.
- Script 05 records how many units fall in a cluster of coincident centroids and how many of those
  clusters mix taxa, beside the thinning result it already exported.
- The repository carries the data. `Dados/`, `08_Dados_Especies/` and `Outputs/` mirror the layout the
  scripts address with `here()`; only the rasters are held back for size.

## Revision of 22 August 2026

Three defects were found by tracing every number in the manuscript back to the cell that produces it,
rather than by matching numbers against the set of all exported values. All three are fixed in the
scripts and documented at the point where they apply.

- Script 02, Section 9. The fractions of `varpart()` were read by position. In vegan 2.6.4 the rows of
  `$part$indfract` are ordered `[a] = X1|X2`, `[b] = X2|X1`, `[c]`, `[d]`, so positional reading
  exchanged the shared and the purely spatial fractions. Two arithmetic identities are now asserted
  before anything is written.
- Script 02, Section 11. The ordination response was the square root of the absolute class area, which
  carries the area of the sampling unit. It is now Hellinger transformed. Term tests were sequential and
  were reported as marginal; marginal tests are now computed and reported.
- Script 03, Section 3. Three input files were accepted in order of preference, of which only one
  reproduces the published results. The required input is now named.

## Citation

Please cite the manuscript once published. Until then, cite this repository and its DOI.

## References

Blanchet, F.G., Legendre, P., Borcard, D. (2008) Forward selection of explanatory variables. *Ecology*
89, 2623-2632.
Brewer, M.J., Butler, A., Cooksley, S.L. (2016) The relative performance of AIC, AICc and BIC in the
presence of unobserved heterogeneity. *Methods in Ecology and Evolution* 7, 679-692.
Brooks, T.M. et al. (2019) Measuring terrestrial Area of Habitat (AOH) and its utility for the IUCN Red
List. *Trends in Ecology & Evolution* 34, 977-986.
Fahrig, L. (2013) Rethinking patch size and isolation effects: the habitat amount hypothesis. *Journal
of Biogeography* 40, 1649-1663.
Gelber, S. et al. (2025) Geometric and demographic effects explain contrasting fragmentation-biodiversity
relationships across scales. *Oikos* 2025, e10778.
Diniz-Filho, J.A.F., Rangel, T.F.L.V.B., Bini, L.M. (2008) Model selection and information theory in
geographical ecology. *Global Ecology and Biogeography* 17, 479-488.
Kronmal, R.A. (1993) Spurious correlation and the fallacy of the ratio standard revisited. *Journal of
the Royal Statistical Society A* 156, 379-392.
Neel, M.C., McGarigal, K., Cushman, S.A. (2004) Behavior of class-level landscape metrics across
gradients of class aggregation and area. *Landscape Ecology* 19, 435-455.
Tredennick, A.T., Hooker, G., Ellner, S.P., Adler, P.B. (2021) A practical guide to selecting models for
exploration, inference, and prediction in ecology. *Ecology* 102, e03336.
Yates, L.A., Aandahl, Z., Richards, S.A., Brook, B.W. (2023) Cross validation for model selection: a
review with examples from ecology. *Ecological Monographs* 93(1), e1557.

## Revisão de 25 de agosto de 2026

A cadeia foi reexecutada pela autora. Dos 66 arquivos CSV em
`Outputs/Manuscrito/`, três mudaram: `TableS28_Geometric_Coupling.csv`, agora
coerente com o script 05 depois da correção de entrada de 24 de agosto, com a
verificação da identidade algébrica passando de 0,029 para 0,0277;
`PartV/TableS_PartV_Identity_and_PartialCorrelations.csv`, só em rótulos; e
`Table_S_PartIII_ObsPred_R2.csv`, arquivo novo. Os outros 63 são idênticos byte a
byte.

As 42 tabelas do material suplementar e as 3 do manuscrito foram conferidas
célula a célula contra os arquivos exportados, com casamento de linha por chave e
exigência de que cada coluna venha de um único campo do CSV em todas as linhas
(`doc/verify_grid2.py`). Trinta e oito passaram automaticamente; as outras quatro
foram conferidas à mão e também estão corretas.

Quatro correções de texto foram aplicadas, nenhuma numérica: duas frases do
manuscrito que contradiziam a tabela que citavam, dois arredondamentos errados na
legenda da Tabela S21, e uma oração sem predicado na legenda da Tabela S28f.
`docs/Registro_de_Correcoes_25.08.md` traz o antes e o depois de cada uma.

Três scripts novos, fora da cadeia inferencial, calculam a probabilidade de
conectividade de Saura e Pascual-Hortal (2007) em três envelopes, como material
de defesa. Ver `docs/TABLE_MAP.md`.

---

## Estado em 29 de agosto de 2026

`Dados/Processados/Data_Raw_WithCoords.csv` foi substituido nesta data. A copia
anterior do repositorio trazia a variavel resposta arredondada para inteiro, com
valores diferentes dos publicados, de modo que uma execucao a partir do
repositorio nao reproduzia os numeros do artigo. O arquivo agora e identico ao da
pasta de trabalho, md5 `881996be7335b4af0c90bcc5c5303e41`.

`R/08_sampling_effort_sensitivity.R` foi executado em R 4.5.2 nesta data. As
Tabelas S42 e S42b do material suplementar vem dessa execucao, e o `session_info`
correspondente esta em `Outputs/Manuscrito/PartVI/session_info_part08.txt` da
pasta de trabalho.

`docs/Registro_de_Correcoes_29.08.md` traz a auditoria numerica completa e o
antes e o depois de cada edicao.
