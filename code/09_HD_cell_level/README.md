This directory contains Visium HD analysis done at the cell level. In other
words, segmentation is performed and counts are aggregated from the original
2um * 2um bins into individual cells for downstream analysis.

## Building a cell-level `SpatialExperiment` and quality control

- `01_bin2cell.*`: Use `bin2cell` to segment Visium HD into cells, aggregate
gene expression into cells, and end up with an `AnnData` in Python that contians
cells by genes
- `exploratory/enact/01_run_enact.*`: Experiment with using the `ENACT` pipeline in place
of the `bin2cell` step
- `02_build_spe_raw.*`: Convert `AnnData` to `SpatialExperiment`, which will be
the format used for the remaining analyses. Don't filter out data or perform
normalization yet
- `03_build_spe_QC.*`: Perform QC, filtering genes and cells. Log normalize
counts
- `07_xenium_genes.*`: Check how thoroughly the Xenium gene panel is expressed
in Visium HD

## Highly variable genes and PCA

- `04_HVG_PCA.*`: Find highly variable genes and add PCA (computed on these
HVGs) to the normalized `SpatialExperiment` in place

## Clustering cells

- `05_banksy_embedding.*`: Compute the `Banksy` embedding for this data, which
is by far the most computationally expensive step in finding spatially informed
cellular clusters with `Banksy`. See `06_banksy_clustering.*` for the rest
- `06_banksy_clustering.*`: Cluster the `Banksy` embedding computed in
`05_banksy_embedding.*` using Leiden clustering at several resolutions

## Spatial registration

See `registration_banksy` for how Banksy clusters found in
`06_banksy_clustering.*` are registered against snRNA-seq data, multi-ome data,
and Visium standard clusters.
