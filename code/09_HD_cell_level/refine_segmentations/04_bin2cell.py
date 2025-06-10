import matplotlib.pyplot as plt
import scanpy as sc
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
out_path = here(
    'processed-data', '09_HD_cell_level', 'probe_fix', 'refine_segmentations',
    'adata', f'{sample_id}_post_QC.h5ad'
)
plot_dir = here(
    'plots', '09_HD_cell_level', 'probe_fix', 'refine_segmentations', 'filtered'
)
mpp = 0.3
expansion_distance = 5
num_random_cells = 3
random_state = 0

os.makedirs(plot_dir, exist_ok=True)

################################################################################
#   Gather and QC segmentations
################################################################################

print(f"{datetime.datetime.now()} | Loading and performing QC on segmentations")

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
#   Plot random regions containing primary segmentations
################################################################################

print(f"{datetime.datetime.now()} | Plotting primary segmentations")

random_cells = (
    adata.obs
        .loc[adata.obs['labels_he'] != 0, :]
        .drop_duplicates(subset = 'labels_he')
        .sample(n = num_random_cells, random_state = random_state)
        ['labels_he']
        .values
)

for i, random_cell in enumerate(random_cells):
    small_adata = adata[adata.obs['labels_he'] == random_cell, :]

    mask = (
        (adata.obs['array_row'] >= small_adata.obs['array_row'].min() - 40) &
        (adata.obs['array_row'] <= small_adata.obs['array_row'].max() + 40) &
        (adata.obs['array_col'] >= small_adata.obs['array_col'].min() - 40) &
        (adata.obs['array_col'] <= small_adata.obs['array_col'].max() + 40) &
        (adata.obs['labels_he'] != 0)
    )
    small_adata = adata[mask, :].copy()
    small_adata.obs['labels_he'] = small_adata.obs['labels_he'].astype(str)

    #   Plot the primary-cell labels
    sc.pl.spatial(
        small_adata, color=[None, "labels_he"],
        img_key=f"{mpp}_mpp_150_buffer", basis="spatial_cropped_150_buffer"
    )
    plt.savefig(
        os.path.join(plot_dir, f'{sample_id}_{i}.png')
    )
    plt.close('all')

################################################################################
#   Aggregate bins into cells
################################################################################

print(f"{datetime.datetime.now()} | Aggregating bins into cells and saving")

adata = b2c.bin_to_cell(
    adata, labels_key="labels_joint",
    spatial_keys=["spatial", "spatial_cropped_150_buffer"]
)
sc.write(out_path, adata)

session_info.show()
