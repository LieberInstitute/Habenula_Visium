import scanpy as sc
import pandas as pd
from pyhere import here
import bin2cell as b2c

sample_id = 'H1-MVPY9BW_A1_8433'
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

#   Read in extracellular bins for this sample
extra_bins = pd.read_csv(extra_bins_path)
extra_bins = extra_bins[extra_bins['sample_id'] == sample_id]

extra_bins['spot_name'] = extra_bins['cell_id'].astype(str) + '_' + extra_bins['sample_id']
extra_bins.set_index('bin_id', inplace=True)

#   Read in habenula-associated bins for this sample
hb_anno = pd.read_csv(hb_anno_path)
hb_anno = hb_anno[hb_anno['sample_id'] == sample_id]

#   Annotate each bin with a cell
adata.obs['spot_name'] = extra_bins['spot_name']
adata.obs.dropna(inplace=True)

#   Grab only bins that belong to habenula
adata.obs['bin_id'] = adata.obs.index
adata.obs.set_index('spot_name', inplace=True)
hb_anno.set_index('spot_name', inplace=True)
adata.obs['hb_anno'] = hb_anno['ManualAnnotation']
adata.obs.dropna(inplace=True)

#   Set index to be bins (which are unique)
adata.obs['spot_name'] = adata.obs.index
adata.obs.set_index('bin_id', inplace=True)

sc.write(adata, out_path)
