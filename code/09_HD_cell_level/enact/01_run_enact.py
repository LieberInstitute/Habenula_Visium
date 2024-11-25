from enact.pipeline import ENACT
from pyhere import here
import os
import session_info

sample_id = 'H1-W369TJK_D1_9090'
bin_to_cell_method = "weighted_by_area"
cell_annotation_method = "celltypist"

out_dir = here('processed-data', '09_HD_cell_level', 'enact', sample_id)
in_dir = here(
    'processed-data', '01_spaceranger', sample_id, 'outs', 'binned_outputs',
    'square_002um'
)
raw_image_path = here('raw-data', 'images', 'vis-hd', f'{sample_id}.tif')
cell_typist_model = here(
    'processed-data', '09_HD_cell_level', 'enact', 'Developing_Human_Brain.pkl'
)

so_hd = ENACT(
    cache_dir = out_dir,
    wsi_path = raw_image_path,
    visiumhd_h5_path = os.path.join(in_dir, 'filtered_feature_bc_matrix.h5'),
    tissue_positions_path = os.path.join(in_dir, 'spatial', 'tissue_positions.parquet'),
    analysis_name = sample_id,
    bin_to_cell_method = bin_to_cell_method,
    cell_annotation_method = cell_annotation_method,
    cell_typist_model = cell_typist_model
)
so_hd.run_enact()

session_info.show()
