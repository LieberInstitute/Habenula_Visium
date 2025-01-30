This directory contains Visium HD analysis done at the cell level. In other
words, segmentation is performed and counts are aggregated from the original
2um * 2um bins into individual cells for downstream analysis.

## Building a cell-level `SpatialExperiment` and quality control

- `01_bin2cell.*`: Use `bin2cell` to segment Visium HD into cells, aggregate
gene expression into cells, and end up with an `AnnData` in Python that contians
cells by genes
- `enact/01_run_enact.*`: Experiment with using the `ENACT` pipeline in place
of the `bin2cell` step
- `02_build_spe_raw.*`: Convert `AnnData` to `SpatialExperiment`, which will be
the format used for the remaining analyses. Don't filter out data or perform
normalization yet
- `03_build_spe_QC.*`: Perform QC, filtering genes and cells. Log normalize
counts
- `04_quality_comparison.*`: Generate QC plots, also comparing QC metrics to
other Visium HD samples (DLPFC, HPC) and Visium standard experiments from the
same brain regions. A more up-to-date version of this script also includes LC
samples and is maintained [here](https://github.com/LieberInstitute/lc_visium_hd/blob/devel/code/04_QC/01_sample_level.R)

## Highly variable genes and PCA

- `05_HVG_PCA.*`: Find highly variable genes and add PCA (computed on these
HVGs) to the normalized `SpatialExperiment` in place
- `06_plot_HVGs.*`: Plot the top 12 HVGs by rank

## Clustering cells

- `07_meringue.*`: Use `MERINGUE` to find spatially informed cellular clusters
- `08_banksy_embedding.*`: Compute the `Banksy` embedding for this data, which
is by far the most computationally expensive step in finding spatially informed
cellular clusters with `Banksy`. See `10_banksy_clustering.*` for the rest
- `10_banksy_clustering.*`: Cluster the `Banksy` embedding computed in
`08_banksy_embedding.*` using k means for many values of k

## Finding spatial expression patterns

- `09_meringue_genes.*`: Explore other functionality from `MERINGUE` (other than
the spatially informed cellular clustering). Mostly focus on summarizing the
spatial patterns from many spatially variable genes

## Ascribing biological meaning to cellular clusters

- `11_annotation_markers.*`: Find one-vs-all marker genes for `Banksy` clusters
at all values of k, using log-fold changes as t stats for a sort of hacked
spatial registration (done in `12_annotation_cor.*`)
- `12_annotation_cor.*`: Using the manually annotated habenula pilot snRNA-seq
cell types as a reference, perform "hacked spatial registration" to attempt to
label `Banksy` clusters with cell types. Plot a heatmap of cell types vs
clusters
- `13_singleR.*`: With a similar goal as the `12_annotation_cor.*` analysis,
directly deconvolve the Visium HD data using `SingleR` and the habenula pilot
snRNA-seq data as a reference
- `14_annotation_comparison.*`: Compare cell-type calls between the annotated
Banksy clusters from  `12_annotation_cor.*` and the deconvolved cells from
`13_singleR.*`
