library(here)
library(tidyverse)
library(spatialLIBD)
library(HDF5Array)
library(sessioninfo)

spe_dir = here(
    'processed-data', '10_HD_bin_level', 'probe_fix', 'spe_norm_filtered'
)
plot_dir = here('plots', '10_HD_bin_level', 'probe_fix', 'marker_genes')
marker_genes = list(
    white_matter = c("MBP", "GFAP", "PLP1", "AQP4"),
    habenula = c("POU4F1", "GPR151", "CHRNB4", "HTR2C"),
    #   originally included "ADARB2", which is no longer expressed (it's
    #   measured in the raw object with 0 counts)
    thalamus = c("LYPD6B", "RORB")
)

dir.create(plot_dir, showWarnings = FALSE)

spe = loadHDF5SummarizedExperiment(spe_dir)
spe$exclude_overlapping = FALSE

#   All markers should be measured
stopifnot(all(unlist(marker_genes) %in% rowData(spe)$gene_name))

#   Convert from gene symbol to Ensembl ID
marker_genes = lapply(
    marker_genes, function(x) rownames(spe)[match(x, rowData(spe)$gene_name)]
)

#   Plot marker combinations for each sample
for (sample_id in unique(spe$sample_id)) {
    for (marker_name in names(marker_genes)) {
        p <- vis_gene(
            spe, sampleid = sample_id, geneid = marker_genes[[marker_name]],
            multi_gene_method = "z_score", is_stitched = TRUE,
            point_size = 10, spatial = TRUE, cap_percentile = 0.995
        ) +
            #   Essentially, make legend bigger while preserving most other
            #   vis_gene defaults for theme
            theme_bw(base_size = 30) +
            theme(
                panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_blank(),
                axis.line = element_blank(),
                axis.text = element_blank(),
                axis.ticks = element_blank(),
                legend.box.spacing = unit(0, "pt")
            )

        png(
            file.path(plot_dir, sprintf("%s_%s.png", marker_name, sample_id)),
            width = 2000, height = 2000
        )
        print(p)
        dev.off()
    }
}

session_info()
