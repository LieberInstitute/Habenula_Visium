This directory contains Visium HD analysis done at the bin level (expression
is a matrix of genes by bins), rather than attempting to form individual cells.
The 8um * 8um bin size is used in most cases, as recommended by 10x Genomics.

## Building a bin-level `SpatialExperiment` and quality control

- `01_build_spe.*`: Read in data as a `SpatialExperiment`, filter out genes and
bins without expression, and perform log normalization
- `02_quality_comparison.*`: Generate QC plots, also comparing QC metrics to
other Visium HD samples (DLPFC, HPC) and Visium standard experiments from the
same brain regions

## Clustering and spatial domain finding

- `03_ficture_transcripts.sh`: Prepare inputs for `FICTURE`
- `04_ficture_run.sh`: Run the full `FICTURE` pipeline for finding subcellular
spatial domains
- `05_hergast.*`: Find spatial domains with `HERGAST`

## Spatially variable genes

- `06_rasterize.*`: Use `SEraster` to lower the resolution of the Visium sample
to approximately Visium-standard resolution. The idea is that SVGs can be
computed with similar accuracy in a computationally reasonable time with lower-
resolution data
- `07_nnSVG.*`: Run `nnSVG` to find spatially variable genes on the lower-
resolution data
- `08_plot_SVGs.*`: Plot the top 12 SVGs by rank

## Cell-cell communication (attempted)

- `09_nest_preprocess.sh`: Run the `NEST` preprocessing step, one of several
steps in a cell-cell communication pipeline. Required prohibitively large
amounts of memory (> 1TB) for this HD data
- `10_nest_run.sh`: Another `NEST` step that was discontinued after the first
`09_nest_preprocess.sh` script failed to run
