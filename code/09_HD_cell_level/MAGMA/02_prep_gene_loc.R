#   Take the GTF used for the multiome data and lift it over to hg19. Everything
#   must be in hg19 for MAGMA (hg38 seems better, but it seems extremely
#   complicated to get the 1000 Genomes European plink files needed for MAGMA in
#   hg38)

library(rtracklayer)
library(here)
library(tidyverse)
library(sessioninfo)

reference_gtf = '/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/10x/refdata-gex-GRCh38-2024-A/genes/genes.gtf.gz'
chain_path = '/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Hb_multiome/processed-data/10_MAGMA/hg38ToHg19.over.chain'
out_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'MAGMA',
    'hg19_gene_loc.tsv'
)

gtf = import(reference_gtf)
gtf = gtf[gtf$type == "gene"]

chain = import.chain(chain_path)
liftOver(gtf, chain) |>
    unlist() |>
    as.data.frame() |>
    as_tibble() |>
    mutate(chr = gsub("chr", "", seqnames)) |>
    #   In the case of duplicate mappings, take the largest interval
    group_by(gene_id) |>
    arrange(desc(end - start)) |>
    slice_head(n = 1) |>
    ungroup() |>
    select(gene_id, chr, start, end) |>
    write_tsv(out_path, col_names = FALSE)

session_info()
