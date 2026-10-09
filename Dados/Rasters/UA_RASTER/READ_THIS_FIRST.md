# The per-unit rasters are not in this repository

Script `00_recompute_class_metrics.R` reads 67 GeoTIFF rasters from this folder, one per sampling unit,
named `ud<N>_<taxon>1.tif` (for example, `ud1_aguariba1.tif`). Together they occupy about 2.6 GB, more
than a Git repository should carry, so they are archived separately and this folder holds only this
note.

Download them from the data repository recorded in `CITATION.cff` and place them here without renaming
them. Script 00 will then run.

No other script reads the rasters. Script 00 recomputes the class-level metrics from them, checks that
the stored values reproduce, and writes `Dados/Processados/Tabela_Metricas_Recomputada.csv`, which is
in the repository. Scripts 01 to 08 start from that table and from the other inputs listed in
`README.md`.
