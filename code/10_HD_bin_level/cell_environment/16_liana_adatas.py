#   For LIANA+, we need cellular and extracellular AnnDatas with more info than
#   what's exported in 03_extracellular_bins.py:
#       - do gene + cell filtering, including QC artifact removal
#       - library size normalization + log transform
#       - add cell-type annotation
#       - add tissue section info
#       - add manual region info (habenula + thalamus)
#       - fix anatomical orientation (spatial coordinates)
#   Produce one dataset-wide AnnData each for cellular + extracellular with
#   this info

import scanpy as sc
import os
from pyhere import here
import session_info
import pandas as pd
import anndata as ad

#   We'll either operate on cellular or extracellular data depending on the
#   array task ID
task_id = int(os.getenv('SLURM_ARRAY_TASK_ID'))
if task_id == 1:
    adata_in_dir = here('processed-data', '09_HD_cell_level', 'new_samples2')
    dataset = 'cellular'
elif task_id == 2:
    adata_in_dir = here(
        'processed-data', '10_HD_bin_level', 'no_secondary', 'cell_environment',
        'adata'
    )
    dataset = 'extracellular'
else:
    raise ValueError(f'Unexpected task_id: {task_id}')

sample_info_path = here('raw-data', 'sample_info', 'hd_basic_info.csv')
tissue_path = here(
    'processed-data', '10_HD_bin_level', 'no_secondary', 'cell_environment',
    'tissue_key_map.csv.gz'
)
hb_anno_path = here(
    'processed-data', '09_HD_cell_level', 'new_samples2',
    'hb_thal_manual_anno.csv.gz'
)
banksy_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'banksy',
    'leiden_res1_8.csv'
)
ct_anno_path = here(
    'processed-data', '09_HD_cell_level', 'no_secondary', 'registration_banksy',
    'cluster_annotation.csv'
)
adata_out_dir = here(
    'processed-data', '10_HD_bin_level', 'no_secondary', 'cell_environment',
    'adata', 'liana_ready'
)
H1_MVPY9BW_A1_8433_artifact = 34326

os.makedirs(adata_out_dir, exist_ok=True)

#   Read in and form one dataset-wide AnnData
adata_list = []
for sample_id in pd.read_csv(sample_info_path)['sample_id']:
    adata = sc.read(os.path.join(adata_in_dir, f'{sample_id}.h5ad'))
    
    #   Only the extracellular object has the 'key' variable so far. It also
    #   needs secondary segmentations dropped
    if dataset == 'cellular':
        adata = adata[adata.obs['labels_joint_source'] == 'primary', :].copy()
        adata.obs['key'] = adata.obs.index + '_' + sample_id
        adata.obs.set_index('key', inplace=True)
    
    adata_list.append(adata)

adata = ad.concat(adata_list, axis=0)

#   Filter out the technical artifact in H1-MVPY9BW_A1_8433
adata = adata[
        ~adata.obs.index.str.contains(r'_8433$') | 
        (adata.obsm['spatial'][:,0] <= H1_MVPY9BW_A1_8433_artifact),
        :
    ].copy()

#   Perform filtering and log normalization similar to that done for the
#   cell-level object in R
sc.pp.filter_genes(adata, min_cells=1)
sc.pp.filter_cells(adata, min_counts=10)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

#   Annotate cells with manual region labels
hb_anno = pd.read_csv(hb_anno_path, index_col='spot_name')
adata.obs['region'] = hb_anno['ManualAnnotation']
adata.obs['region'] = adata.obs['region'].fillna('other').astype('category')

#   Add cluster and cell type labels (drop the tiny fraction of cells missing
#   these)
banksy_df = pd.read_csv(banksy_path, index_col='key')
ct_anno_df = pd.read_csv(ct_anno_path, index_col='cluster')

adata.obs['banksy_cluster'] = banksy_df['banksy']
adata.obs['cell_type'] = adata.obs['banksy_cluster'].map(
    ct_anno_df['fine_cell_type']
)

missing_prop = adata.obs['cell_type'].isna().mean()
print(f'Percentage of cells missing a cell-type label: {missing_prop:.3%}')
adata = adata[~adata.obs['cell_type'].isna(), :].copy()

tissue_df = pd.read_csv(tissue_path, index_col='key')
adata.obs['tissue_section'] = tissue_df['tissue_section']
adata.obsm['spatial'] = tissue_df.loc[
    adata.obs_names, ['pxl_col_in_fullres', 'pxl_row_in_fullres']
].to_numpy()
missing_prop = adata.obs['tissue_section'].isna().mean()
print(f'Percentage of cells missing a tissue section label: {missing_prop:.3%}')
if missing_prop > 0:
    adata = adata[~adata.obs['tissue_section'].isna(), :].copy()

sc.write(os.path.join(adata_out_dir, f'{dataset}.h5ad'), adata)

session_info.show()
