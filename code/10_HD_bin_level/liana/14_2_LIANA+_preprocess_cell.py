# This script matches single-cell RNA-seq data to spatial transcriptomics data using nearest neighbor search.
# It reads in spatial data and cell annotations, finds the nearest spatial spot for each cell, and merges the data.
# Finally, it saves the merged AnnData object for further analysis.
# Here I used cell level data to match with cell level data

import scanpy as sc
import os
import pandas as pd
from pyhere import here
import bin2cell as b2c
import numpy as np
from sklearn.neighbors import KDTree
from scipy.spatial import cKDTree
import matplotlib.pyplot as plt

extra_bins_path = here(
    'processed-data', '10_HD_bin_level', 'no_secondary', 'cell_environment',
    'extracellular_bins.csv.gz'
)

adata = sc.read_h5ad("/dcs04/lieber/lcolladotor/Habenula_R01_LIBD4270/Habenula_Visium/processed-data/10_HD_bin_level/no_secondary/cell_environment/adata/liana_ready/cellular.h5ad")
adata = adata[adata.obs['region'] == 'habenula', :].copy()

out_path = here(
    'processed-data', '10_HD_bin_level', "no_secondary",'liana','adata'
)
out_path.mkdir(parents=True, exist_ok=True)

# file_path = os.path.join(out_path, f"adata_cellular_withcelltype.h5ad")
# adata = sc.read(file_path)

gtf_path = "/dcs04/lieber/lcolladotor/annotationFiles_LIBD001/10x/refdata-gex-GRCh38-2020-A/genes/genes.gtf"

gtf = pd.read_csv(
    gtf_path,
    sep='\t',
    comment='#',
    header=None,
    names=["chrom", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"]
)

gtf_genes = gtf[gtf["feature"] == "gene"]

def parse_attributes(attr_str):
    attrs = {}
    for item in attr_str.strip().split(';'):
        if item.strip():
            key, value = item.strip().split(' ', 1)
            attrs[key] = value.strip('"')
    return attrs

attrs_parsed = gtf_genes["attribute"].apply(parse_attributes)
attr_df = pd.DataFrame(attrs_parsed.tolist())

gene_map = (attr_df[["gene_id","gene_name"]]
            .drop_duplicates()
            .assign(gene_id=lambda df: df["gene_id"].str.replace(r"\.\d+$", "", regex=True)))

adata.var["var_names_backup"] = adata.var_names.astype(str)

if "gene_id" not in adata.var.columns:
    adata.var["gene_id"] = adata.var_names.to_series()

adata.var["gene_id"] = adata.var["gene_id"].astype(str).str.replace(r"\.\d+$", "", regex=True)

if "gene_name" in adata.var.columns:
    adata.var = adata.var.drop(columns=["gene_name"])

adata.var = adata.var.merge(gene_map, on="gene_id", how="left")

name = adata.var["gene_name"].astype("string")
gid  = adata.var["gene_id"].astype("string")

missing = name.isna() | (name.str.lower().isin(["nan", "none", ""]))
name_filled = name.where(~missing, gid)

adata.var_names = pd.Index(name_filled.astype(str))
adata.var_names_make_unique(join="-")

n_total = adata.var.shape[0]
n_unmapped = missing.sum()
print(f"[INFO] total genes: {n_total}, newly unmapped-to-symbol (filled by gene_id): {int(n_unmapped)}")
print(adata.var_names[:20])

if "gene_name" in adata.var.columns:
    adata.var["gene_name_orig"] = adata.var["gene_name"]
    adata.var = adata.var.drop(columns=["gene_name"])

adata.write(os.path.join(out_path, f"adata_cellular_withcelltype.h5ad"))

print("✅ Saved:", os.path.join(out_path, "adata_cellular_withcelltype.h5ad"))