import matplotlib.pyplot as plt
import scanpy as sc
import os
import pandas as pd
from pyhere import here
import session_info
import bin2cell as b2c
import datetime

sample_id = 'H1-W369TJK_D1_9090'
stardist_dir = here('processed-data', '09_HD_cell_level', 'stardist')
final_out_path = here('processed-data', '09_HD_cell_level', f'{sample_id}.h5ad')
pre_out_path = here(
    'processed-data', '09_HD_cell_level', f'{sample_id}_pre_bin2cell.h5ad'
)
sr_dir = here(
    'processed-data', '01_spaceranger', sample_id, 'outs', 'binned_outputs',
    'square_002um'
)
sr_spatial_dir = here(
    'processed-data', '01_spaceranger', sample_id, 'outs', 'spatial'
)
plot_dir = here('plots', '09_HD_cell_level')
raw_image_path = here('images', 'vis-hd', f'{sample_id}.tif')
mpp = 0.3

os.makedirs(stardist_dir, exist_ok=True)
os.makedirs(plot_dir, exist_ok=True)

################################################################################
#   Build and preprocess AnnData
################################################################################

print(f"{datetime.datetime.now()} | Building and preprocessing AnnData")

#   Read in spaceranger outputs into an AnnData
adata = b2c.read_visium(
    sr_dir,
    source_image_path = raw_image_path,
    spaceranger_image_path = sr_spatial_dir
)

#   Use Ensembl IDs for var_names
adata.var_names = adata.var['gene_ids']

#   Require bins with nonzero counts and genes present in at least 3 bins
sc.pp.filter_genes(adata, min_cells=3)
sc.pp.filter_cells(adata, min_counts=1)

#   Create a scaled H&E image attached to the object (and for segmentation with
#   stardist)
b2c.scaled_he_image(
    adata,
    mpp = mpp,
    save_path = os.path.join(stardist_dir, f'he_{sample_id}.tiff')
)

#   Normalize counts to account for "striping" effect
b2c.destripe(adata)

################################################################################
#   Perform nuclear-based ("primary") segmentation
################################################################################

print(f"{datetime.datetime.now()} | Performing nuclear-based ('primary') segmentation")

#   Segment nuclei on H&E image
b2c.stardist(
    image_path=os.path.join(
        stardist_dir, f'he_{sample_id}.tiff'
    ),
    labels_npz_path=os.path.join(
        stardist_dir, f'he_{sample_id}.npz'
    ),
    stardist_model="2D_versatile_he", 
    prob_thresh=0.01
)

#   Add segmentations to object
b2c.insert_labels(
    adata, 
    labels_npz_path=os.path.join(
        stardist_dir, f'he_{sample_id}.npz'
    ), 
    basis="spatial", 
    spatial_key="spatial_cropped",
    mpp=mpp, 
    labels_key="labels_he"
)

#   Expand labels to attempt to capture cells and not nuclei
b2c.expand_labels(
    adata, 
    labels_key='labels_he', 
    expanded_labels_key="labels_he_expanded"
)

################################################################################
#   Perform gene-expression-based ('secondary') segmentation
################################################################################

print(f"{datetime.datetime.now()} | Performing gene-expression-based ('secondary') segmentation")

#   Create an image from gene counts
b2c.grid_image(
    adata,
    "n_counts_adjusted",
    mpp=mpp,
    sigma=5,
    save_path=os.path.join(
        stardist_dir, f'gex_{sample_id}.tiff'
    )
)

#   Segment cells on the gene-count image
b2c.stardist(
    image_path=os.path.join(
        stardist_dir, f'gex_{sample_id}.tiff'
    ), 
    labels_npz_path = os.path.join(
        stardist_dir, f'gex_{sample_id}.npz'
    ), 
    stardist_model="2D_versatile_fluo", 
    prob_thresh=0.05, 
    nms_thresh=0.5
)

#   Add segmentations to object
b2c.insert_labels(
    adata, 
    labels_npz_path = os.path.join(
        stardist_dir, f'gex_{sample_id}.npz'
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

#-------------------------------------------------------------------------------
#   Plot primary and secondary cells
#-------------------------------------------------------------------------------

#   Region for plots
mask = (
    (adata.obs['array_row'] >= 1000) & 
    (adata.obs['array_row'] <= 1050) & 
    (adata.obs['array_col'] >= 1000) & 
    (adata.obs['array_col'] <= 1050)
)

#   Plot union of cell labels
bdata = adata[mask]
bdata = bdata[bdata.obs['labels_joint'] > 0]
bdata.obs['labels_joint'] = bdata.obs['labels_joint'].astype(str)
sc.pl.spatial(
    bdata, color=[None, "labels_joint_source", "labels_joint"],
    img_key=f"{mpp}_mpp", basis="spatial_cropped"
)
plt.savefig(
    os.path.join(plot_dir, f'{sample_id}_cells.png')
)
plt.close('all')

#   Keep a copy of the AnnData before aggregation (to enable interactive
#   plotting later, for example)
sc.write(pre_out_path, adata)

################################################################################
#   Aggregate bins into cells
################################################################################

print(f"{datetime.datetime.now()} | Aggregating bins into cells")

adata = b2c.bin_to_cell(
    adata, labels_key="labels_joint",
    spatial_keys=["spatial", "spatial_cropped"]
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
    bdata, color="bin_count", img_key=f"{mpp}_mpp", basis="spatial_cropped"
)
plt.savefig(
    os.path.join(plot_dir, f'{sample_id}_cells_aggregated.png')
)
plt.close('all')

sc.write(final_out_path, adata)
session_info.show()
