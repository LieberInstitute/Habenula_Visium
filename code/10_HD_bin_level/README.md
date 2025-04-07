This directory contains Visium HD analysis done at the bin level (expression
is a matrix of genes by bins), rather than attempting to form individual cells.
The 8um * 8um bin size is used in most cases, as recommended by 10x Genomics.

## Building a bin-level `SpatialExperiment` and quality control

- `01_build_spe.*`: Read in data as a `SpatialExperiment`, filter out genes and
bins without expression, and perform log normalization
- `02_QC.*`: Identify which bins to drop based on quality metrics

## Spatially variable genes

- `03_rasterize.*`: Use `SEraster` to lower the resolution of the Visium sample
to approximately Visium-standard resolution. The idea is that SVGs can be
computed with similar accuracy in a computationally reasonable time with lower-
resolution data
- `04_nnSVG.*`: Run `nnSVG` to find spatially variable genes on the lower-
resolution data
- `05_gather_variable_genes.*`: Gather nnSVG results from each sample to compute
dataset-wide SVGs. Export, then plot top SVGs and HVGs

## Clustering

See `ficture_harmony` for clustering with FICTURE.

## Other

- `06_plot_markers.*`: Plot habenula, thalamus, and white-matter markers for
each sample
