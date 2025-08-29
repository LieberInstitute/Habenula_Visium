#Change the data to AnnDate format

#--------------------------
# conda activate giotto_env
# pip install --upgrade pip
# pip install anndata scanpy

library(reticulate)
use_python("/jhpce/shared/libd/core/visium_hd/1.0/hd_env/bin/python", required = TRUE)
GiottoClass::set_giotto_python_path("/jhpce/shared/libd/core/visium_hd/1.0/hd_env/bin/python")
reticulate::py_config()

export RETICULATE_PYTHON=/users/cliu3/.conda/envs/giotto_env/bin/python
export RETICULATE_PYTHON=/jhpce/shared/libd/core/visium_hd/1.0/hd_env/bin/python
#--------------------------
library(HDF5Array)
library(SingleCellExperiment) # or SummarizedExperiment, if not SCE
library(spatialLIBD)  # for loadHDF5SummarizedExperiment, if custom, use as needed
library(Giotto)
library(zellkonverter)

spe_dir <- "/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/09_HD_cell_level/probe_fix/spe_norm_filtered"
spe <- loadHDF5SummarizedExperiment(spe_dir)

sample_ids <- unique(colData(spe)$sample_id)

for (sid in sample_ids) {
    message("Processing sample: ", sid)
    
    # Subset the spe object for this sample
    spe_sub <- spe[, colData(spe)$sample_id == sid]
    
    # Convert to Giotto object
    g_obj <- spatialExperimentToGiotto(spe_sub)
     
    # Convert Giotto object to AnnData (zellkonverter expects SCE, so use Giotto -> SCE -> AnnData if needed)
    # If giottoToAnnData is available:
    adata <- giottoToAnnData(g_obj)
    
    # Save AnnData file
    anndata_file <- paste0("/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/LIANA/cell_level_", sid, ".h5ad")
    writeH5AD(adata, file = anndata_file)
    
    message("Saved AnnData for: ", sid)
}

