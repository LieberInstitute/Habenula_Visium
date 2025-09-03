import pandas as pd
import scanpy as sc
import os
from pyhere import here
import session_info
import bin2cell as b2c
import anndata as ad
import matplotlib.pyplot as plt

import sys
sys.path.append(str(here('code', '10_HD_bin_level', 'cell_environment')))
import extracellular_bins_functions as ebf

ad_in_paths = here(
    'processed-data', '09_HD_cell_level', 'probe_fix',
    '{}.h5ad'
)
ad_pre_in_paths = here(
    'processed-data', '09_HD_cell_level', 'probe_fix',
    '{}_pre_bin2cell.h5ad'
)
hb_anno_path = here(
    'processed-data', '09_HD_cell_level', 'probe_fix', 'crawdad', 'region',
    'habenula_shiny_annotations.csv.gz'
)
out_dir = here(
    'processed-data', '09_HD_cell_level', 'probe_fix', 'LIANA', 'adata'
)
plot_dir = here('plots', '09_HD_cell_level', 'probe_fix', 'LIANA')
min_bins_per_cell = 4
expansion_distance = 5

os.makedirs(plot_dir, exist_ok=True)
os.makedirs(out_dir, exist_ok=True)

################################################################################
#   Ordinary cell-level data
################################################################################

sample_id_path = here('raw-data', 'sample_info', 'hd_sample_list.txt')
with open(sample_id_path, 'r') as f:
    all_samples = f.read().splitlines()[:3]

adata_list = []
for sample_id in all_samples:
    adata = sc.read(str(ad_in_paths).format(sample_id))
    adata.obs['sample_id'] = sample_id
    adata.obs['key'] = adata.obs.index + '_' + adata.obs['sample_id']
    adata.obs.index = adata.obs['key']
    adata_list.append(adata)

adata = ad.concat(adata_list, axis=0)

#   Perform filtering and log normalization similar to that done for the
#   cell-level object in R
sc.pp.filter_genes(adata, min_cells=1)
sc.pp.filter_cells(adata, min_counts=10)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

#   Filter out the artifact in H1-MVPY9BW_A1_8433
adata = adata[
    (adata.obs['sample_id'] != 'H1-MVPY9BW_A1_8433') |
    (adata.obsm['spatial'][:,0] <= 34223),
    :
]

#   Annotate cells as habenula or not
hb_anno = pd.read_csv(hb_anno_path, index_col='spot_name')
adata.obs['region'] = hb_anno['ManualAnnotation']
adata.obs['region'] = adata.obs['region'].fillna('other').astype('category')

sc.write(os.path.join(out_dir, 'cellular.h5ad'), adata)

################################################################################
#   Extracellular data
################################################################################

adata_list = []
for sample_id in all_samples:
    adata = sc.read(str(ad_pre_in_paths).format(sample_id))

    adata = ebf.find_microenvironment(
        adata, expansion_distance = expansion_distance
    )

    #   Add 'cell_component' column, needed for 'ebf.drop_bad_secondary_cells'
    adata.obs['cell_component'] = 'Unlabeled'
    adata.obs.loc[adata.obs['microenvironment_secondary'] != 0, 'cell_component'] = 'Sec. Extracellular'
    adata.obs.loc[adata.obs['labels_gex'] != 0, 'cell_component'] = 'Sec. Cell Body'
    adata.obs.loc[adata.obs['microenvironment_primary'] != 0, 'cell_component'] = 'Prim. Extracellular'
    adata.obs.loc[adata.obs['labels_he_expanded'] != 0, 'cell_component'] = 'Prim. Cell Body'
    adata.obs.loc[adata.obs['labels_he'] != 0, 'cell_component'] = 'Prim. Nucleus'

    #   Surprisingly, bin2cell passively drops many secondary cells by giving all of
    #   their bins (either this or none of the bins) a 'none' value in
    #   'labels_joint_source'. Just drop such bins now (cell body only)
    adata = adata[
        (adata.obs['labels_joint_source'] == 'secondary') |
        (adata.obs['cell_component'] != 'Sec. Cell Body')
    ]

    adata = ebf.drop_bad_secondary_cells(
        adata, min_bins_per_cell = min_bins_per_cell
    )

    #   Since some originally secondary cells are dropped, numbering is thrown
    #   off relative to the cells we already have annotated. Create a map from
    #   microenvironment cell labels to the original cell labels
    cell_map = (
        adata.obs
            #   Since 'primary' comes before 'secondary', this makes it so
            #   originally secondary cells aren't later labeled with a primary
            #   label that "takes over" through expansion
            .sort_values(by = 'microenvironment_joint_source', ascending=False)
            .loc[
                adata.obs['labels_joint'] != 0,
                ['microenvironment_joint', 'labels_joint']
            ]
            .drop_duplicates(subset = 'microenvironment_joint', keep = 'first')
            .set_index('microenvironment_joint')
    )
    cell_map.index = cell_map.index.astype(str)

    #   Drop intracellular bins
    adata = adata[
        adata.obs['cell_component'].isin(
            ['Prim. Extracellular', 'Sec. Extracellular']
        )
    ]

    adata = b2c.bin_to_cell(
        adata, labels_key="microenvironment_joint",
        spatial_keys=["spatial", "spatial_cropped_150_buffer"]
    )

    #   Bring in original cell labels
    adata.obs['labels_joint'] = cell_map['labels_joint']
    assert all(~adata.obs['labels_joint'].isna())

    #   Add sample ID and label cells with original labels
    adata.obs['sample_id'] = sample_id
    adata.obs['key'] = adata.obs['labels_joint'].astype(str) + '_' + adata.obs['sample_id']
    adata.obs.set_index('key', inplace=True)

    adata_list.append(adata)

adata = ad.concat(adata_list, axis=0)

#   Perform filtering and log normalization similar to that done for the
#   cell-level object in R
sc.pp.filter_genes(adata, min_cells=1)
sc.pp.filter_cells(adata, min_counts=10)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

#   Filter out the artifact in H1-MVPY9BW_A1_8433
adata = adata[
    (adata.obs['sample_id'] != 'H1-MVPY9BW_A1_8433') |
    (adata.obsm['spatial'][:,0] <= 34223),
    :
]

adata.obs['region'] = hb_anno['ManualAnnotation']
adata.obs['region'] = adata.obs['region'].fillna('other').astype('category')

#   For one sample, visually validate habenula region labels
small_adata = adata[adata.obs['sample_id'] == 'H1-MVPY9BW_A1_8433']
sc.pl.spatial(
    small_adata, color=["region"], basis="spatial_cropped_150_buffer",
    spot_size=200.0
)
plt.savefig(os.path.join(plot_dir, 'H1-MVPY9BW_A1_8433_hb_cells.png'))
plt.close('all')

sc.write(os.path.join(out_dir, 'extracellular.h5ad'), adata)

session_info.show()
