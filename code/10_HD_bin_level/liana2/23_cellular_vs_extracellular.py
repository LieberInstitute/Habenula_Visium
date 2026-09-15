# Compare extracellular and cellular heatmap raw or normalized scores in a scatter plot to identify pairs of cells types that communicate differently
# read the data from heatmap

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
from pyhere import here

plot_dir1= here(
    'plots', '10_HD_bin_level', 'no_secondary', 'liana2'
)
plot_dir2= here(
    'plots', '10_HD_bin_level', 'no_secondary', 'liana2', 'extracellular'
)
os.makedirs(plot_dir1, exist_ok=True)
os.makedirs(plot_dir2, exist_ok=True)   

bandwidth = 5000.0
cellular_heatmap_path = os.path.join(plot_dir1,f"source_target_sum_mean_heatmap_data_{bandwidth}.csv")
extracellular_heatmap_path = os.path.join(plot_dir2,f"source_target_sum_mean_heatmap_data_{bandwidth}.csv")
cellular_heatmap = pd.read_csv(cellular_heatmap_path, index_col=0)
extracellular_heatmap = pd.read_csv(extracellular_heatmap_path, index_col=0)

