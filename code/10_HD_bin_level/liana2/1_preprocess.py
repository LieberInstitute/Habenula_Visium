# add the cell type annotation to the cell-level data

import scanpy as sc
import os
import pandas as pd
from pyhere import here
import bin2cell as b2c
import numpy as np

out_path = here(
    'processed-data', '10_HD_bin_level', "new_samples2",'liana', 'input_habenula',
    'cellular_annotated.h5ad'
)  
anno_data_path = here(
    'processed-data', '10_HD_bin_level', 'no_secondary','cell_environment','adata',
    'liana_ready','cellular.h5ad'
)
anno_data = sc.read_h5ad(anno_data_path)

# read in the cell type info
cell_type_info_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'registration_banksy','cluster_annotation.csv')

cell_type_info = pd.read_csv(cell_type_info_path)

# --- 1. combine the data ---
anno_data.obs["banksy_cluster"] = anno_data.obs["banksy_cluster"].astype("Int64")
cell_type_info["cluster"] = cell_type_info["cluster"].astype("Int64")

anno_data.obs = (
    anno_data.obs
    .reset_index()
    .merge(
        cell_type_info,
        left_on="banksy_cluster",
        right_on="cluster",
        how="left"
    )
    .set_index("key")
)

# --- 1. filter out excit.thal ---
anno_data = anno_data[anno_data.obs["fine_cell_type"]!= "Excit.Thal"]

# --- 2. Save the filtered AnnData ---
sc.write(out_path, anno_data)