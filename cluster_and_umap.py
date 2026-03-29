# -*- coding: utf-8 -*-
"""

"""

#%%

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

# %%

data_df = pd.read_csv('MC_SC_17_s03L.csv')
cell_id = data_df["Cell ID"].values
data_df = data_df.iloc[:,7:23]

data_df['slc17a6'] = data_df['slc17a6'].astype(np.float64)
data_df['slc32a1'] = data_df['slc32a1'].astype(np.float64)
data_df['rph3a'] = data_df['rph3a'].astype(np.float64)
data_df['slc5a7'] = data_df['slc5a7'].astype(np.float64)

adata = ad.AnnData(data_df)
adata.obs_names = cell_id

# %%

sc.pp.log1p(adata)
sc.pp.scale(adata)
sc.tl.pca(adata)

sc.pp.neighbors(adata)
sc.tl.umap(adata)

# %%

sc.pl.umap(adata,color='phox2b')

