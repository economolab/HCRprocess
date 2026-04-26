# -*- coding: utf-8 -*-
"""

"""

import anndata as ad
import numpy as np
import pandas as pd

# %%

filepaths = ['MC_SC_17_s03L.csv', 'MC_SC_17_s01R.csv']

usecols = ['Cell ID', 'slc17a6', 'slc32a1', 'slc5a7', 'phox2b', 
           'ralyl', 'tenm2', 'ebf3', 'pcp4', 'tshz2', 'alcam', 'celf2', 
           'meis2', 'rph3a', 'robo1', 'syt1', 'zfhx3']

#%%

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

def build_cell_anndata(cell_df):
    
    genes = list(cell_df.columns)
    genes = [gene.capitalize() for gene in genes]
    cell_df.columns = genes
    
    for gene in genes:
        cell_df[gene] = cell_df[gene].astype(np.uint16)
        
    adata = ad.AnnData(cell_df)
        
    return adata

# %%

cell_dfs = []

for filepath in filepaths:
    cell_dfs.append(load_cell_df(filepath, usecols, index_col='Cell ID'))
    
cell_df = concat_cell_dfs(cell_dfs)
cell_adata = build_cell_anndata(cell_df)

cell_adata.write_h5ad("MapMyCellsInput.h5ad")

# %% load example MapMyCells input file

import scanpy as sc

example_file = sc.read_h5ad(r"MapMyCells/wholemousebrain_ccn20230722_example_10kcells_550genes.h5ad")