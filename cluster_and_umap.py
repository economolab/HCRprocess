# -*- coding: utf-8 -*-
"""

"""

#%%

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

# %%

filepaths = ['MC_SC_17_s03L.csv', 'MC_SC_17_s01R.csv']

usecols = ['Cell ID', 'egfp', 'dtom', 'slc17a6', 'slc32a1', 'slc5a7', 'phox2b', 
           'ralyl', 'tenm2', 'ebf3', 'pcp4', 'tshz2', 'alcam', 'celf2', 
           'meis2', 'rph3a', 'robo1', 'syt1', 'zfhx3', 'Area']

obs_anno = ['egfp', 'dtom']

# %%

def load_cell_df(filepath, usecols, index_col='Cell ID'):
    
    cell_df = pd.read_csv(filepath, index_col=index_col, usecols=usecols)
    
    for col in cell_df.columns:
        if not isinstance(cell_df[col].dtype, np.float64):
            cell_df[col] = cell_df[col].astype(np.float64)
    
    return cell_df

def concat_cell_dfs(cell_dfs):
    
    columns = [list(cell_df.columns) for cell_df in cell_dfs]
    if not all(cols == columns[0] for cols in columns):
        print("Columns do not match across all dataframes, concatenation failed.")
        return
    
    cell_df = pd.concat(cell_dfs)
    
    return cell_df

def build_cell_anndata(cell_df, obs_anno):
    
    obs_anno = dict.fromkeys(obs_anno)
    
    for key in obs_anno.keys():
        obs_anno[key] = cell_df[key].values
        cell_df.drop(columns=key, inplace=True)
        
    adata = ad.AnnData(cell_df)
    
    adata.layers['raw'] = np.copy(adata.X)
    
    for key in obs_anno.keys():
        adata.obs[key] = obs_anno[key]
        
    return adata

# %%

cell_dfs = []

for filepath in filepaths:
    cell_dfs.append(load_cell_df(filepath, usecols))
    
cell_df = concat_cell_dfs(cell_dfs)

adata = build_cell_anndata(cell_df, obs_anno)
    
# %%

sc.pp.log1p(adata)
sc.pp.scale(adata)
sc.tl.pca(adata)

sc.pp.neighbors(adata)
sc.tl.umap(adata, min_dist=0.5)
sc.tl.leiden(adata, flavor="igraph", n_iterations=-1)

# %%

sc.pl.umap(adata, color='ralyl')

# %%

import matplotlib.pyplot as plt

raw_exp = np.copy(adata.layers['raw'])
area = raw_exp[:,-1]
area = np.expand_dims(area, axis=1)
raw_exp = np.delete(raw_exp, -1, axis=1)
raw_exp = raw_exp / area

R = np.corrcoef(np.transpose(raw_exp))

fig, ax = plt.subplots()
im = ax.imshow(R)

genes = adata.var_names
genes = genes[0:-1]
ax.set_xticks(range(len(genes)), labels=genes,
              rotation=45, ha="right", rotation_mode="anchor")
ax.set_yticks(range(len(genes)), labels=genes)

for i in range(len(genes)):
    for j in range(len(genes)):
        text = ax.text(j, i, np.round(R[i, j], decimals=1),
                       ha="center", va="center", color="w")
        
# %%

import os
from ABC_toolbox import cell_funcs, gene_funcs

local_data_dir = r'C:\Users\jpv88\Documents\GenePicker9001_data'
abc_genes = np.load(r'C:\Users\jpv88\Documents\GitHub\GenePicker9001\ABC_toolbox\util_files\gene_scRNAseq.npy', 
                    allow_pickle=True)
abc_genes = [gene.lower() for gene in abc_genes]
abc_genes_idx = [np.where(gene == np.array(abc_genes))[0][0] for gene in genes]

meta = pd.read_csv(os.path.join(local_data_dir, "antIRN-PARN-scRNAseq-meta.csv"), low_memory=False)
exp = np.load(os.path.join(local_data_dir, "antIRN-PARN-scRNAseq-raw.npy"))
freqs = pd.read_pickle(os.path.join(local_data_dir, "antIRN-PARN-MERFISH-freqs.pkl"))
exp = exp[:,abc_genes_idx]

exp_super, meta_super = cell_funcs.boot_super(exp, meta, k=3)

exp_super = gene_funcs.normalize_counts_to_median(exp_super)

# bootstrap distribution from scRNAseq that matches MERFISH frequencies
exp_boot, meta_boot = cell_funcs.bootstrap_scRNAseq(meta_super, exp_super, freqs, 
                                                    n=10000)

# %%

R = np.corrcoef(np.transpose(exp_boot))

fig, ax = plt.subplots()
im = ax.imshow(R)

genes = adata.var_names
genes = genes[0:-1]
ax.set_xticks(range(len(genes)), labels=genes,
              rotation=45, ha="right", rotation_mode="anchor")
ax.set_yticks(range(len(genes)), labels=genes)

for i in range(len(genes)):
    for j in range(len(genes)):
        text = ax.text(j, i, np.round(R[i, j], decimals=1),
                       ha="center", va="center", color="w")

# %%

for gene in data_df.columns:

    sc.pl.umap(adata, color=gene)
    
# %%

sc.pl.umap(adata, color='egfp')

# %%

sc.pl.pca(adata, color='slc17a6')

