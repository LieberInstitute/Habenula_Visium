#   For each HD sample, try different probability thresholds for the secondary
#   segmentations and visually examine

import matplotlib.pyplot as plt
import scanpy as sc
import os
from pyhere import here
import session_info
import bin2cell as b2c
import datetime
import pandas as pd

#   Use the task ID to loop over sample and probability threshold combinations
task_id = int(os.getenv('SLURM_ARRAY_TASK_ID')) - 1
sample_index = task_id % 5
prob_index = task_id // 5

sample_info_path = here('raw-data', 'sample_info', 'hd_basic_info.csv')
sample_info = pd.read_csv(sample_info_path)
sample_id = sample_info.iloc[sample_index]['sample_id']
spaceranger_dir = sample_info.iloc[sample_index]['spaceranger_dir']
batch_num = sample_info.iloc[sample_index]['batch_num']
prob_thres = round(0.2 + 0.1 * prob_index, 1)

stardist_new_dir = here(
    'processed-data', '09_HD_cell_level', 'new_samples2',
    'refine_segmentations', 'stardist'
)
stardist_old_dir = here(
    'processed-data', '09_HD_cell_level', 'new_samples2', 'stardist'
)
pre_out_path = here(
    'processed-data', '09_HD_cell_level', 'new_samples2',
    f'{sample_id}_pre_bin2cell.h5ad'
)
out_path = here(
    'processed-data', '09_HD_cell_level', 'new_samples2',
    'refine_segmentations', 'adata',
    f'{sample_id}_{str(prob_thres).replace(".", "_")}.h5ad'
)
plot_dir = here(
    'plots', '09_HD_cell_level', 'new_samples2', 'refine_segmentations',
    str(prob_thres).replace('.', '_')
)
mpp = 0.3

os.makedirs(stardist_new_dir, exist_ok=True)
os.makedirs(plot_dir, exist_ok=True)
os.makedirs(out_path.parent, exist_ok=True)

################################################################################
#   Build and preprocess AnnData
################################################################################

print(f"{datetime.datetime.now()} | Loading AnnData")

adata = sc.read_h5ad(pre_out_path)

################################################################################
#   Perform gene-expression-based ('secondary') segmentation
################################################################################

print(f"{datetime.datetime.now()} | Performing gene-expression-based ('secondary') segmentation")

#   Segment cells on the gene-count image
b2c.stardist(
    image_path=os.path.join(
        stardist_old_dir, f'gex_{sample_id}.tiff'
    ), 
    labels_npz_path = os.path.join(
        stardist_new_dir, f'gex_{sample_id}.npz'
    ), 
    stardist_model="2D_versatile_fluo", 
    prob_thresh=prob_thres, 
    nms_thresh=0.5
)

#   Add segmentations to object
b2c.insert_labels(
    adata, 
    labels_npz_path = os.path.join(
        stardist_new_dir, f'gex_{sample_id}.npz'
    ), 
    basis="array", 
    mpp=mpp, 
    labels_key="labels_gex"
)

#   Take the union of cell labels from both segmentation methods
b2c.salvage_secondary_labels(
    adata, 
    primary_label="labels_he_expanded", 
    secondary_label="labels_gex", 
    labels_key="labels_joint"
)

sc.write(out_path, adata)

################################################################################
#   Plot primary and secondary cells
################################################################################

#   Plot 2 different subregions to get a representative idea
for i in range(2):
    #   Region for plots
    mask = (
        (adata.obs['array_row'] >= 1000 + 500 * i) & 
        (adata.obs['array_row'] <= 1050 + 500 * i) & 
        (adata.obs['array_col'] >= 1000 + 500 * i) & 
        (adata.obs['array_col'] <= 1050 + 500 * i)
    )

    #   If the region has no cells, try to iterate over other regions until
    #   cells are found
    offset = 150
    while ((adata[mask].obs['labels_he'] == 0).sum() < 10) and (offset < 1000):
        mask = (
            (adata.obs['array_row'] >= 1000 + 500 * i + offset) & 
            (adata.obs['array_row'] <= 1050 + 500 * i + offset) & 
            (adata.obs['array_col'] >= 1000 + 500 * i + offset) & 
            (adata.obs['array_col'] <= 1050 + 500 * i + offset)
        )
        offset += 150
    assert (adata[mask].obs['labels_he'] == 0).sum() >= 10, "Failed to find a region with cells for plotting"
    
    print(f'Image {i+1} bounds: ({1000 + 500 * i + offset}, {1050 + 500 * i + offset})')

    #   Plot primary segmentations
    crop = b2c.get_crop(
        adata[mask], basis="spatial", spatial_key="spatial_cropped_150_buffer",
        mpp=mpp
    )
    rendered = b2c.view_labels(
        image_path = os.path.join(
            stardist_old_dir, f'he_{sample_id}.tiff'
        ),
        labels_npz_path = os.path.join(
            stardist_old_dir, f'he_{sample_id}.npz'
        ),  
        crop = crop
    )
    plt.imshow(rendered)
    plt.savefig(os.path.join(plot_dir, f'{sample_id}_primary_{i+1}.png'))
    plt.close('all')

    #   Plot secondary segmentations
    crop = b2c.get_crop(adata[mask], basis="array", mpp=mpp)
    rendered = b2c.view_labels(
        image_path = os.path.join(
            stardist_old_dir, f'gex_{sample_id}.tiff'
        ),
        labels_npz_path = os.path.join(
            stardist_new_dir, f'gex_{sample_id}.npz'
        ),  
        crop = crop,
        stardist_normalize = True
    )
    plt.imshow(rendered)
    plt.savefig(os.path.join(plot_dir, f'{sample_id}_secondary_{i+1}.png'))
    plt.close('all')

################################################################################
#   Check secondary proportion
################################################################################

adata = b2c.bin_to_cell(
    adata, labels_key="labels_joint",
    spatial_keys=["spatial", "spatial_cropped_150_buffer"]
)

print(f"Secondary proportion: {(adata.obs['labels_joint_source'] == 'secondary').mean():.2f}")

session_info.show()
