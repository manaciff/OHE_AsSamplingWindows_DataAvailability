# Where each table and figure comes from

The files written by the scripts are named after the part of the chain that
produces them. The manuscript and the supplementary material number their tables
in reading order. The two numbering schemes do not coincide, so this file maps
one onto the other. Every number reported in the two documents can be traced to
a row of a file listed here.

Paths are relative to `outputs/`, which is `Outputs/Manuscrito/` in the working
folder of the project.

## Main text

| Item | Produced by | File |
|---|---|---|
| Figure 1, study area map | ArcGIS Pro, outside the R chain | not in the repository; see Supplementary Section 5 |
| Figure 2, composition ridgeline | `01` | `PartI/Figure2_Ridgeline_Composition.tiff` |
| Figure 3, configuration ridgeline | `01` | `PartI/Figure3_Ridgeline_Configuration.tiff` |
| Figure 4, NMDS | `01` | `PartI/Figure1_NMDS_Composition_Genus.tiff` |
| Figure 5, RDA biplot | `02` | `PartII/Figure_RDA_Biplot.tiff` |
| Figure 6, coefficients of the reported model | `05` | `PartV/Figure_PV_02_Reported_Model_Coefficients.tiff` |
| Figure 7, genus-level heatmap | `03` | `Figure_PIII_08_Effect_Heatmap.tiff` |
| Figure 8, conservation value index | ArcGIS Pro, outside the R chain | not in the repository |
| Table 1, global PERMANOVA | `01` | `PartI/TableS4_PERMANOVA_Global.csv` |
| Table 2, marginal tests of the RDA | `02` | `PartII/TableII_4_RDA_Anova_Marginal.csv` |
| Table 3, reported model of envelope extent | `05` | `PartV/TableS_PartV_Uncoupled_Model.csv` |

## Supplementary material

| Table | Produced by | File |
|---|---|---|
| S1, S2 | written by hand | no producing script |
| S3, S4 | `01` | `PartI/Table3_Medians_Composition_PLAND.csv`, `PartI/Table3b_Medians_MinorClasses.csv`, `PartI/Table1_Composition_By_Species.csv`, `PartI/Table1b_MinorClasses_And_LPI_By_Species.csv` |
| S5, S6 | `01` | `PartI/Table4_Medians_Configuration.csv`, `PartI/Table2_Configuration_By_Species.csv` |
| S7, envfit | `01` | `PartI/TableS2_Envfit_Configuration.csv` |
| S8, pairwise PERMANOVA by genus | `01` | `PartI/TableS5_PERMANOVA_Pairwise_Genus.csv` |
| S9, VIF of the RDA pool | `02` | `PartII/TableII_1_VIF_RDA.csv` |
| S10, VIF of the partitioning pool | `02` | `PartII/TableII_1_VIF_VP.csv` |
| S11, betadisper | `01` | `PartI/TableS3_Betadisper.csv` |
| S12, S13, S14, pairwise PERMANOVA | `01` | `PartI/TableS5_PERMANOVA_Pairwise_Species.csv`, `..._Diet.csv`, `..._Locomotion.csv` |
| S15, descriptive statistics | `03` | `Table_S_PartIII_DescriptiveStatistics.csv` |
| S16, comparison of the six specifications | `03` | `Table_S_PartIII_ModelComparison.csv` |
| S17, residual diagnostics | `03` | `Table_S_PartIII_Diagnostics.csv` |
| S18, VIF of the a priori model | `03` | `Table_S_PartIII_VIF.csv` |
| S18a, full a priori model | `03` | `Table_S_PartIII_ConfirmatoryModel.csv` |
| S18b, model averaging | `03` | `Table_S_PartIII_Averaging_FullVsConditional.csv` |
| S19, selection table | `03` | `Table_S_PartIII_Dredge_FullTable.csv` |
| S20, variable importance | `03` | `Table_S_PartIII_VariableImportance.csv` |
| S21, stability check | `03` | `Table_S_PartIII_StabilityCheck.csv` |
| S22, FRAC against patch area | `03` | `Table_S_PartIII_FRAC_Area_Correlation.csv` |
| S23, matrix sensitivity | `03` | `Table_S_PartIII_Sensitivity_MatrixModel.csv` |
| S24, genus-level regressions | `03` | `Table_S_PartIII_EffectByGenus.csv` |
| S25, genus-level Spearman | `03` | `Table_S_PartIII_MetricAssociation_ByGenus.csv` |
| S26, mixed-model coupling controls | `04` | `TableS28_Geometric_Coupling.csv` |
| S27, sample sizes | `01`, `02` | `PartI/TableS_Sample_Sizes_PartI.csv`, `PartII/TableII_0_Sample_Sizes.csv` |
| S28a, six specifications | `05` | `PartV/TableS_PartV_Specifications.csv` |
| S28b, permutation null | `05` | `PartV/TableS_PartV_NullSimulation.csv` |
| S28c, AICc under the null | `05` | `PartV/TableS_PartV_AICc_Under_Null.csv` |
| S28d, decision table | `05` | `PartV/TableS_PartV_Decision.csv` |
| S28e, algebraic coupling | `05` | `PartV/TableS_PartV_Identity_and_PartialCorrelations.csv` |
| S28f, gain from the interaction | `05` | `PartV/TableS_PartV_Interaction_Gain.csv` |
| S29, genus as moderator | `03` | `Table_S_PartIII_SlopeByGenus_Moderator.csv` |
| S30, proximity against amount | `05` | `PartV/TableS_PartV_Proximity_Amount_Correlations.csv` |
| S31, proximity conditioned | `05` | `PartV/TableS_PartV_Proximity_Conditioned.csv` |
| S32, agreement of the recomputation | `00` | `Part0/TableS_Recompute_Agreement.csv` |
| S33, repaired cells | `00` | `Part0/TableS_Repaired_Cells.csv` |
| S34, spatial sensitivity | `03` | `Table_S_PartIII_SpatialSensitivity.csv` |
| S35, criteria and out-of-sample error | `05` | `PartV/TableS_PartV_CV_And_BIC.csv` |
| S36, S37, conservation value index | ArcGIS Pro, outside the R chain | not in the repository |
| S38, partition without PD | `02` | `PartII/TableII_3b_VP_Fractions_Without_PD.csv` |
| S39, effect of the ordination transformation | `02` | `PartII/TableII_4b_RDA_Transformation_Comparison.csv`, `PartII/TableII_5_Response_Size_Dependence.csv` |
| S40, farming matrix in an uncoupled specification | `05` | `PartV/TableS_PartV_Matrix_Uncoupled.csv` |
| S41, thinning test for spatial coincidence | `05` | `PartV/TableS_PartV_Spatial_Thinning.csv` |
| S42, sampling effort bounds | `08` | `PartVI/TableS42_Sampling_Effort.csv` |
| S42b, sampling effort robustness | `08` | `PartVI/TableS42b_Sampling_Effort_Robustness.csv` |
| Figure S1, correlation matrix | `02` | `PartII/FigS1_Collinearity_Pearson.tiff` |
| Figure S2, response transformation | `03` | `Figure_PIII_01_Response_Transformation.tiff` |
| Figure S3, variable importance | `03` | `Figure_PIII_03_Variable_Importance.tiff` |
| Figure S4, partial effects | `03` | `Figure_PIII_04_Partial_Effects.tiff` |
| Figure S5, observed against predicted | `03` | `Figure_PIII_05_Observed_vs_Predicted.tiff`; the two values quoted in the caption are in `Table_S_PartIII_ObsPred_R2.csv` |
| Figure S6, residual diagnostics | `03` | `Figure_PIII_06_Residual_Diagnostics.tiff` |
| Figure S7, genus-level Spearman heatmap | `03` | `Figure_PIII_07_Metric_Association_Heatmap.tiff` |
| Figure S8, variation partitioning | `02` | `PartII/Figure_VP_Venn.tiff` |
| Figure S9, permutation null | `05` | `PartV/Figure_PV_01_Null_Distribution.tiff` |
| Figure S10, coefficients of the a priori set | `03` | `Figure_PIII_02_Coefficient_ForestPlot.tiff` |

## Files that no longer correspond to anything reported

`PartII/TableII_4_RDA_Anova_Terms.csv` was renamed to
`TableII_4_RDA_Anova_Terms_Sequential.csv` on 22 August 2026, when the
manuscript moved from sequential to marginal term tests. Delete the old file
before archiving so that no reader takes it for the reported table.

`Figure_PIII_02_Coefficient_ForestPlot.tiff` is Figure S10 of the supplement
since 24 August 2026. It was written by script `03` and reproduced nowhere, which
is why the entry above used to say it corresponded to nothing. Figure 6 of the
manuscript shows the coefficients of the reported model and comes from script
`05`; Figure S10 shows the six terms of the a priori set, with the three the
response normalises drawn hollow and unmarked.

## Which model the figures of Part III show

Section 8 of script `03` chooses the retained specification among the models
within two AICc units, first by the smallest number of parameters and then by the
smaller cross-validated error. On these data that is the Gamma model with a
logarithmic link, so `global_model` is a `glm` and every figure downstream is
drawn from it. Figures S5 and S6 are therefore panels of the Gamma model, not of
the mixed one. Their captions said otherwise until 23 August 2026.

## A demonstração de conectividade não está arquivada aqui

Um cálculo exploratório da probabilidade de conectividade de Saura e Pascual-Hortal (2007) e da
importância de cada mancha, dPC, dentro de alguns envelopes foi preparado em 25 de agosto de 2026 como
material de defesa. Não aparece no manuscrito nem no material suplementar, e **nem os scripts nem as
figuras estão arquivados neste repositório**. O registro fica aqui pelo argumento, que toca a Seção
2.5.1 do manuscrito.

PC é normalizado pela área total da paisagem, que neste desenho é a área do envelope e portanto a
variável resposta, de modo que PC não pode ser comparado entre unidades amostrais. dPC é a queda
relativa em PC quando uma mancha é removida, e a área da paisagem ao quadrado cancela entre numerador e
denominador. O ordenamento das manchas dentro de uma mesma unidade é, por isso, livre do acoplamento
geométrico; um ordenamento de unidades inteiras por PC não seria.
