import matplotlib.pyplot as plt
import scanpy as sc
import pandas as pd
import numpy as np
import os
from pyhere import here
import session_info
import bin2cell as b2c
import datetime

sample_id_path = here('raw-data', 'sample_info', 'hd_sample_list.txt')
with open(sample_id_path, 'r') as f:
    all_samples = f.read().splitlines()
sample_id = all_samples[int(os.getenv('SLURM_ARRAY_TASK_ID')) - 1]

orig_path = here(
    'processed-data', '09_HD_cell_level', 'probe_fix',
    f'{sample_id}_pre_bin2cell.h5ad'
)
plot_dir = here(
    'plots', '09_HD_cell_level', 'probe_fix', 'refine_segmentations', 'filtered'
)
mpp = 0.3
expansion_distance = 5

os.makedirs(plot_dir, exist_ok=True)

################################################################################
#   Something
################################################################################

adata = sc.read(orig_path)

#   Label microenvironment around primary segmentations (nuclei)
b2c.expand_labels(
    adata, 
    labels_key='labels_he_expanded', 
    expanded_labels_key="microenvironment_primary",
    max_bin_distance = expansion_distance
)

#   Drop primary cells that have no surrounding bins (to control overly dense
#   segmentations)
mask = ~ adata.obs['labels_he_expanded'].isin(
    adata.obs.loc[
        adata.obs['labels_he_expanded'] == 0, 'microenvironment_primary'
    ]
)
adata.obs.loc[mask, 'labels_he_expanded'] = 0

b2c.salvage_secondary_labels(
    adata, 
    primary_label="labels_he_expanded", 
    secondary_label="labels_gex", 
    labels_key="labels_joint"
)

################################################################################
#   Aggregate bins into cells
################################################################################

print(f"{datetime.datetime.now()} | Aggregating bins into cells")

adata = b2c.bin_to_cell(
    adata, labels_key="labels_joint",
    spatial_keys=["spatial", "spatial_cropped_150_buffer"]
)

cell_mask = (
    (adata.obs['array_row'] >= 1450) & 
    (adata.obs['array_row'] <= 1550) & 
    (adata.obs['array_col'] >= 250) & 
    (adata.obs['array_col'] <= 450)
)

#   Plot counts within cells after aggregation of bins
bdata = adata[cell_mask]
sc.pl.spatial(
    bdata, color="bin_count", img_key=f"{mpp}_mpp_150_buffer",
    basis="spatial_cropped_150_buffer"
)
plt.savefig(
    os.path.join(plot_dir, f'{sample_id}_cells_aggregated.png')
)
plt.close('all')


sc.write(final_out_path, adata)
session_info.show()
