library(here)
library(tidyverse)
library(SpatialExperiment)
library(spatialLIBD)
library(HDF5Array)
library(sessioninfo)
library(SingleR)

#   Determine variables related to snRNA-seq label resolution
label_names = c("final_Annotations", "final_Annotations_broad")
res_names = c("fine", "broad")
this_label = label_names[as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))]
this_res = res_names[as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))]

sce_path = "/dcs04/lieber/lcolladotor/pilotHb_LIBD001/Roche_Habenula/processed-data/04_snRNA-seq/sce_objects/sce_final.Rdata"
spe_dir = here('processed-data', '10_HD_bin_level', 'spe_norm')
out_path = here(
    'processed-data', '10_HD_bin_level', 'singler', sprintf('%s.csv', this_res)
)

dir.create(dirname(out_path), showWarnings = FALSE)

#   Load unlabeled cell-level SPE and the snRNA-seq reference
load(sce_path, verbose = TRUE)
spe = loadHDF5SummarizedExperiment(spe_dir)

#   Subset each object to the genes in common
shared_genes = intersect(rownames(sce_final), rownames(spe))
spe = spe[shared_genes,]
sce_final = sce_final[shared_genes,]

sce_final$final_Annotations_broad = sub(
    "\\.[0-9]+", "", sce_final$final_Annotations
)

#   Apply annotations and save to CSV
SingleR(
        test = spe, ref = sce_final, labels = sce_final[[this_label]], 
        de.method = "wilcox"
    ) |>
    as_tibble() |>
    write_csv(out_path)

session_info()