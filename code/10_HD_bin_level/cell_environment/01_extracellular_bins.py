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
# sample_id = all_samples[int(os.getenv('SLURM_ARRAY_TASK_ID')) - 1]
sample_id = all_samples[0]

pre_out_path = here(
    'processed-data', '09_HD_cell_level', 'probe_fix',
    f'{sample_id}_pre_bin2cell.h5ad'
)
plot_dir = here('plots', '10_HD_bin_level', 'probe_fix', 'cell_environment')
mpp = 0.3
random_state = 0
num_random_cells = 5

os.makedirs(plot_dir, exist_ok=True)

adata = sc.read(pre_out_path)

#   Label microenvironment around primary segmentations (nuclei)
b2c.expand_labels(
    adata, 
    labels_key='labels_he_expanded', 
    expanded_labels_key="microenvironment_primary",
    max_bin_distance = 3
)

#   Sanity checks: cell labels should be preserved when expanding
#   microenvironment. Microenvironment and cell labels should be one to one
mask = adata.obs['labels_he_expanded'] != 0
assert all(adata.obs['labels_he_expanded'][mask] == adata.obs['microenvironment_primary'][mask]), "Primary microenvironment labels don't always match their original cells"
assert len(adata.obs['labels_he_expanded'].unique()) == len(adata.obs['microenvironment_primary'].unique()), "Primary microenvironment labels should be one to one with cell labels"

#   Label microenvironment around secondary segmentations
b2c.expand_labels(
    adata, 
    labels_key='labels_gex', 
    expanded_labels_key="microenvironment_secondary",
    max_bin_distance = 5
)

#   Sanity checks: cell labels should be preserved when expanding
#   microenvironment. Microenvironment and cell labels should be one to one
mask = adata.obs['labels_gex'] != 0
assert all(adata.obs['labels_gex'][mask] == adata.obs['microenvironment_secondary'][mask]), "Secondary microenvironment labels don't always match their original cells"
assert len(adata.obs['labels_gex'].unique()) == len(adata.obs['microenvironment_secondary'].unique()), "Secondary microenvironment labels should be one to one with cell labels"

b2c.salvage_secondary_labels(
    adata, 
    primary_label="microenvironment_primary", 
    secondary_label="microenvironment_secondary", 
    labels_key="microenvironment_joint"
)

#   Ensure I understand how salvaging works in combination with the
#   microenvironment expansion-- expansion may overwrite previously secondary
#   bins, but never previously primary ones
assert all(adata.obs['microenvironment_joint_source'][adata.obs['labels_joint_source'] == 'primary'] == 'primary')
assert not all(adata.obs['microenvironment_joint_source'][adata.obs['labels_joint_source'] == 'secondary'] == 'secondary')

primary_df = adata.obs[['microenvironment_primary']][
    (adata.obs['microenvironment_joint_source'] == 'primary') &
    (adata.obs['labels_he_expanded'] == 0) &
    (adata.obs['microenvironment_joint'] != 0) # not actually sure why this isn't redundant
]
primary_df = (
    primary_df
        .reset_index()
        .rename(
            {'index': 'bin_id', 'microenvironment_primary': 'cell_id'}, axis = 1
        )
        .assign(sample_id = sample_id)
)

secondary_df = adata.obs[['microenvironment_secondary']][
    (adata.obs['microenvironment_joint_source'] == 'secondary') &
    (adata.obs['labels_gex'] == 0) &
    (adata.obs['microenvironment_joint'] != 0) # not actually sure why this isn't redundant
]
secondary_df = (
    secondary_df
        .reset_index()
        .rename(
            {'index': 'bin_id', 'microenvironment_primary': 'cell_id'}, axis = 1
        )
        .assign(sample_id = sample_id)
)

#   Add 'cell_component' column for informative coloring of plots
adata.obs['cell_component'] = 'Unlabeled'
adata.obs.loc[adata.obs['microenvironment_secondary'] != 0, 'cell_component'] = 'Sec. Extracellular'
adata.obs.loc[adata.obs['labels_gex'] != 0, 'cell_component'] = 'Sec. Cell Body'
adata.obs.loc[adata.obs['microenvironment_primary'] != 0, 'cell_component'] = 'Extracellular'
adata.obs.loc[adata.obs['labels_he_expanded'] != 0, 'cell_component'] = 'Cell Body'
adata.obs.loc[adata.obs['labels_he'] != 0, 'cell_component'] = 'Nucleus'

random_cells = primary_df.sample(
    n = num_random_cells, random_state = random_state
)['cell_id'].values

for i in range(num_random_cells):
    small_adata = adata[adata.obs['microenvironment_primary'] == random_cells[i], :]
    small_adata = adata[
        (adata.obs['array_row'] >= small_adata.obs['array_row'].min() - 40) &
        (adata.obs['array_row'] <= small_adata.obs['array_row'].max() + 40) &
        (adata.obs['array_col'] >= small_adata.obs['array_col'].min() - 40) &
        (adata.obs['array_col'] <= small_adata.obs['array_col'].max() + 40),
        :
    ]
    
    sc.pl.spatial(
        small_adata, color=[None, "cell_component"],
        img_key=f"{mpp}_mpp_150_buffer", basis="spatial_cropped_150_buffer"
    )
    plt.savefig(
        os.path.join(plot_dir, f'{sample_id}_primary_cell{i+1}.png')
    )
    plt.close('all')
