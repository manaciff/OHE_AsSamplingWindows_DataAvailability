# The per-unit rasters are not in this repository

Script `00_recompute_class_metrics.R` reads 67 GeoTIFF rasters from this folder,
one per sampling unit, named `ud<N>_<taxon>1.tif`. Together they occupy about
2.6 GB, which is beyond what a Git repository should carry, so they are archived
separately and this folder is left empty on purpose.

Download them from the Zenodo deposit recorded in `CITATION.cff` and place them
here, keeping the file names unchanged. Script 00 will then run.

You do not need them for anything else. Scripts 01 to 05 read only
`Dados/Processados/`, and every table and figure in the manuscript and in the
supplementary material is produced by those five scripts. Script 00 exists to
recompute the class-level metrics from the rasters and to verify that the stored
values reproduce; `Outputs/Manuscrito/Part0/TableS_Recompute_Agreement.csv` in
this repository is the result of that verification.
