import matplotlib.pyplot as plt
import scanpy as sc
import os
from pyhere import here
import session_info
import bin2cell as b2c

prob_tag = "recommended"

sample_id_path = here('raw-data', 'sample_info', 'hd_sample_list.txt')
with open(sample_id_path, 'r') as f:
    all_samples = f.read().splitlines()
sample_id = all_samples[int(os.getenv('SLURM_ARRAY_TASK_ID')) - 1]

orig_path = here(
    'processed-data', '09_HD_cell_level', 'probe_fix',
    f'{sample_id}_pre_bin2cell.h5ad'
)
plot_dir = here(
    'plots', '09_HD_cell_level', 'probe_fix', 'refine_segmentations', prob_tag
)
stardist_dir = here(
    'processed-data', '09_HD_cell_level', 'probe_fix', 'stardist'
)
mpp = 0.3
num_random_cells = 5
random_state = 0

os.makedirs(plot_dir, exist_ok=True)

################################################################################
#   Plot random regions containing cells
################################################################################

adata = sc.read(orig_path)

random_cells = (
    adata.obs
        .loc[adata.obs['labels_he'] != 0, :]
        .drop_duplicates(subset = 'labels_he')
        .sample(n = num_random_cells, random_state = random_state)
        ['labels_he']
        .values
)

for i, random_cell in enumerate(random_cells):
    small_adata = adata[adata.obs['labels_he'] == random_cell, :].copy()

    mask = (
        (adata.obs['array_row'] >= small_adata.obs['array_row'].min() - 40) &
        (adata.obs['array_row'] <= small_adata.obs['array_row'].max() + 40) &
        (adata.obs['array_col'] >= small_adata.obs['array_col'].min() - 40) &
        (adata.obs['array_col'] <= small_adata.obs['array_col'].max() + 40) &
        adata.obs['labels_he'] != 0
    )
    small_adata = adata[mask, :]
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

    crop = b2c.get_crop(
        adata[mask], basis="spatial", spatial_key="spatial_cropped_150_buffer",
        mpp=mpp
    )
    rendered = b2c.view_labels(
        image_path = os.path.join(
            stardist_dir, f'he_{sample_id}.tiff'
        ),
        labels_npz_path = os.path.join(
            stardist_dir, f'he_{sample_id}.npz'
        ),  
        crop = crop
    )
    #   Plot the primary segmentations
    plt.imshow(rendered)
    plt.savefig(os.path.join(plot_dir, f'{sample_id}_seg{i}.png'))
    plt.close('all')

session_info.show()
