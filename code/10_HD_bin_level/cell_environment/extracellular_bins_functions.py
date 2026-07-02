#   Functions for 01_extracellular_bins.py

import bin2cell as b2c
import matplotlib.pyplot as plt
import scanpy as sc
import os
import pandas as pd

#   Find microenvironment around primary segmentations, expanding around the
#   cell bodies by [expansion_distance] bins. Modify the AnnData in place
def find_microenvironment(adata, expansion_distance):
    #   Label microenvironment around primary segmentations (nuclei)
    b2c.expand_labels(
        adata, 
        labels_key='labels_he_expanded', 
        expanded_labels_key="microenvironment_primary",
        max_bin_distance = expansion_distance
    )

    #   Sanity checks: cell labels should be preserved when expanding
    #   microenvironment. Microenvironment and cell labels should be one to one
    mask = adata.obs['labels_he_expanded'] != 0
    assert all(adata.obs['labels_he_expanded'][mask] == adata.obs['microenvironment_primary'][mask]), "Primary microenvironment labels don't always match their original cells"
    assert len(adata.obs['labels_he_expanded'].unique()) == len(adata.obs['microenvironment_primary'].unique()), "Primary microenvironment labels should be one to one with cell labels"

#   Visualize cell segmentations and surrounding microenvironment. Return a 
#   DataFrame containing just extracellular bins, ready for export
def export_and_plot(adata, plot_dir, sample_id, mpp, random_state = 0):
    #---------------------------------------------------------------------------
    #   Form DataFrame of extracellular bins for export
    #---------------------------------------------------------------------------

    extracellular_df = (
        adata
            .obs[adata.obs['cell_component'] == 'Prim. Extracellular']
            .reset_index(names = 'bin_id')
            .copy()
    )
    extracellular_df['cell_key'] = extracellular_df['microenvironment_primary'].astype(str) + '_' + sample_id

    #---------------------------------------------------------------------------
    #   Visualize cell segmentations and surrounding microenvironment
    #---------------------------------------------------------------------------

    random_cell = (
        extracellular_df
            .sample(n = 1, random_state = random_state)
            ['microenvironment_primary']
            .values[0]
    )
    small_adata = adata[adata.obs['microenvironment_primary'] == random_cell, :]

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
        os.path.join(plot_dir, f'{sample_id}_random_cells.pdf')
    )
    plt.close('all')

    extracellular_df = extracellular_df[['bin_id', 'cell_key']]

    return extracellular_df
