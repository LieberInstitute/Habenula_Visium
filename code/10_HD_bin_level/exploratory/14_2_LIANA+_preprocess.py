import scanpy as sc
import os
import pandas as pd
from pyhere import here
import bin2cell as b2c
import numpy as np
from sklearn.neighbors import KDTree
from scipy.spatial import cKDTree
import matplotlib.pyplot as plt

sample_id_path = here('raw-data', 'sample_info', 'hd_sample_list.txt')
with open(sample_id_path, 'r') as f:
    all_samples = f.read().splitlines()
    sample_id = all_samples[int(os.getenv('SLURM_ARRAY_TASK_ID')) - 1]

sr_dir = here(
    'processed-data', '01_spaceranger', 'probe_fix', sample_id, 'outs',
    'binned_outputs', 'square_002um'
)
sr_spatial_dir = here(
    'processed-data', '01_spaceranger', 'probe_fix', sample_id, 'outs',
    'spatial'
)
raw_image_path = here('raw-data', 'images', 'vis-hd', f'{sample_id}.tif')

hb_anno_path = here(
    'processed-data', '09_HD_cell_level', 'probe_fix', 'crawdad', 'region',
    'habenula_shiny_annotations.csv.gz'
)
extra_bins_path = here(
    'processed-data', '10_HD_bin_level', 'probe_fix', 'cell_environment',
    'extracellular_bins.csv.gz'
)
adata = b2c.read_visium(
    sr_dir,
    count_file = 'filtered_feature_bc_matrix.h5',
    source_image_path = raw_image_path,
    spaceranger_image_path = sr_spatial_dir
)
cell_df_path = here(
    'processed-data', '10_HD_bin_level', 'LIANA'
)
out_path = here(
    'processed-data', '10_HD_bin_level', 'LIANA'
)
cell_df = pd.read_csv(os.path.join(cell_df_path, "cell_df.csv"))
cell_df = cell_df[cell_df["sample_id"].isin([sample_id])].copy()
# match arry_row and array_col in adata.obs with cell_df
obs_xy   = adata.obs[["array_row","array_col"]].astype(float).to_numpy()
cell_xy  = cell_df[["array_row","array_col"]].to_numpy()

# create KDTree and find the nearest spot for each cell
tree = cKDTree(obs_xy)
dist, idx = tree.query(cell_xy, k=1, p=2)  # p=2: Euclidean distance (p=1: Manhattan distance)

# write the coordinates of the nearest spot
cell_df["matched_array_row"] = obs_xy[idx, 0]
cell_df["matched_array_col"] = obs_xy[idx, 1]
cell_df["min_dist"] = dist
cell_df["match_index"] = adata.obs.index[idx]

# show the result
print(cell_df[["key", "cell_type", "array_row","array_col","matched_array_row", "matched_array_col", "min_dist"]].head())

dupe_spots = (
    cell_df
    .groupby(["array_row","array_col","matched_array_row", "matched_array_col"])
    .size()
    .reset_index(name="count")
    .query("count > 1")
)
print(dupe_spots.head())

tree = KDTree(obs_xy)
dist, idx = tree.query(cell_xy, k=1)
cell_df["match_index"] = adata.obs.index[idx.flatten()]

print(cell_df[["key", "cell_type","array_row","array_col", "matched_array_row", "matched_array_col", "min_dist"]].head())
thresh = 0.5
cell_df_filt = cell_df.loc[dist.flatten() <= thresh].copy()

mask = pd.Series(True, index=cell_df.index) 
cell_matched = cell_df.loc[mask, ["match_index", "cell_type", "key", "sample_id"]]

cell_df = cell_df.sort_values("min_dist")

# delete matched spot duplicates, only keep the nearest
cell_df_nodupe = cell_df.drop_duplicates(
    subset=["matched_array_row", "matched_array_col"], keep="first"
).copy()

print(f"after removing duplicates, cell number: {cell_df_nodupe.shape[0]}")

# merge

# construct a DataFrame, with (matched_array_row, matched_array_col) as index
cell_df_nodupe_idx = cell_df_nodupe.set_index(
    ["matched_array_row", "matched_array_col"]
)[["cell_type", "key", "sample_id"]]
cell_df_nodupe_idx.index.set_names(["array_row", "array_col"], inplace=True)

# adata.obs also needs the same MultiIndex
obs_with_index = adata.obs.copy()
obs_with_index = obs_with_index.assign(
    array_row=obs_with_index["array_row"].astype(int),
    array_col=obs_with_index["array_col"].astype(int)
)
obs_with_index.index = pd.MultiIndex.from_arrays(
    [obs_with_index["array_row"], obs_with_index["array_col"]],
    names=["array_row", "array_col"]
)

# combine
adata.obs = obs_with_index.join(cell_df_nodupe_idx, how="left")
adata = adata[~adata.obs["cell_type"].isna()].copy()

adata.write(os.path.join(out_path, f"adata/adata_{adata.obs['sample_id'].unique()[0]}.h5ad"))
