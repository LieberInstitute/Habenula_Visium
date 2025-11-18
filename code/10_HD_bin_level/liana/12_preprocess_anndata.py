# Choose habenula bins from extracellular bins and annotate with habenula region
import scanpy as sc
import os
import pandas as pd
from pyhere import here
import bin2cell as b2c
import numpy as np

sample_id_path = here('raw-data', 'sample_info', 'hd_sample_list.txt')
with open(sample_id_path, 'r') as f:
    all_samples = f.read().splitlines()
sample_id = all_samples[int(os.getenv('SLURM_ARRAY_TASK_ID')) - 1]

out_path = here(
    'processed-data', '10_HD_bin_level', 'nest','input_habenula',
    f'{sample_id}.h5ad'
)

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

# #   Read in extracellular bins for this sample
# extra_bins = pd.read_csv(extra_bins_path)
# extra_bins = extra_bins[extra_bins['sample_id'] == sample_id]

# extra_bins['spot_name'] = extra_bins['cell_id'].astype(str) + '_' + extra_bins['sample_id']
# extra_bins.set_index('bin_id', inplace=True)

# #   Read in habenula-associated bins for this sample
# hb_anno = pd.read_csv(hb_anno_path)
# hb_anno = hb_anno[hb_anno['sample_id'] == sample_id]

# #   Annotate each bin with a cell
# adata.obs['spot_name'] = extra_bins['spot_name']
# adata.obs = adata.obs.join(extra_bins['spot_name'])
# mask = ~adata.obs['spot_name'].isna()
# adata = adata[mask, :].copy()

# #   Grab only bins that belong to habenula
# adata.obs['bin_id'] = adata.obs.index
# hb_anno = hb_anno.set_index('spot_name')
# adata.obs = adata.obs.set_index('spot_name').join(hb_anno['ManualAnnotation'])
# mask = ~adata.obs['ManualAnnotation'].isna()
# adata = adata[mask, :].copy()

# #   Set index to be bins (which are unique)
# adata.obs['spot_name'] = adata.obs.index
# adata.obs = adata.obs.set_index('bin_id')

# assert adata.obs.shape[0] == adata.X.shape[0], (
#     f"obs has {adata.obs.shape[0]} rows but X has {adata.X.shape[0]}"
# )

# sc.write(out_path,adata)

# --- 1. Read in extracellular bins for this sample ---
extra_bins = pd.read_csv(extra_bins_path)
extra_bins = extra_bins[extra_bins['sample_id'] == sample_id].copy()
extra_bins['spot_name'] = extra_bins['cell_id'].astype(str) + '_' + extra_bins['sample_id']
extra_bins.set_index('bin_id', inplace=True)

# --- 2. Read in habenula-associated bins for this sample ---
hb_anno = pd.read_csv(hb_anno_path)
hb_anno = hb_anno[hb_anno['sample_id'] == sample_id].copy()

# --- 3. Annotate each bin with a cell via map ---
adata.obs['spot_name'] = adata.obs_names.map(
    lambda b: extra_bins.loc[b, 'spot_name'] if b in extra_bins.index else np.nan
)

# Drop bins without spot_name (synchronously updates .obs and .X)
adata = adata[adata.obs['spot_name'].notna(), :].copy()

# --- 4. Grab only bins that belong to habenula via map ---
# Build mapping dictionary
hb_dict = hb_anno.set_index('spot_name')['ManualAnnotation'].to_dict()
adata.obs['hb_anno'] = adata.obs['spot_name'].map(hb_dict)

# Drop non-habenula bins
adata = adata[adata.obs['hb_anno'].notna(), :].copy()

# --- 5. Optional: keep spot_name column and check alignment ---
assert adata.obs.shape[0] == adata.X.shape[0], (
    f"obs has {adata.obs.shape[0]} rows but X has {adata.X.shape[0]}"
)

# --- 6. Save the filtered AnnData ---
sc.write(out_path, adata)