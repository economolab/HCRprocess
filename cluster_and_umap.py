# -*- coding: utf-8 -*-
"""

"""

#%%

import anndata as ad
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

mpl.rcParams['image.composite_image'] = False
plt.rcParams['svg.fonttype'] = 'none'

# %%

# filepaths = ['MC_SC_17_s03L.csv', 'MC_SC_17_s01R.csv']

# usecols = ['Cell ID', 'egfp', 'dtom', 'slc17a6', 'slc32a1', 'slc5a7', 'phox2b', 
#            'ralyl', 'tenm2', 'ebf3', 'pcp4', 'tshz2', 'alcam', 'celf2', 
#            'meis2', 'rph3a', 'robo1', 'syt1', 'zfhx3', 'Area']

# obs_anno = ['egfp', 'dtom']

filepaths = ['MC_SC_17_s03L_new_spot.csv']

usecols = ['Cell ID', 'egfp_b', 'dtom_b', 'slc17a6_b', 'slc32a1_b', 'slc5a7_b', 'phox2b', 
           'ralyl', 'tenm2', 'ebf3', 'pcp4', 'tshz2', 'alcam', 'celf2', 
           'meis2', 'rph3a', 'robo1', 'syt1', 'zfhx3', 'Area']

obs_anno = ['egfp_b', 'dtom_b']

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

def build_cell_anndata(cell_df, obs_anno, norm_by_area=False):
    
    obs_anno = dict.fromkeys(obs_anno)
    
    for key in obs_anno.keys():
        obs_anno[key] = cell_df[key].values
        cell_df.drop(columns=key, inplace=True)
    
    if norm_by_area == True:
        
        one_pix = 0.325 * 0.325
        
        for col in cell_df.columns:
            cell_df[col] = cell_df[col] / (cell_df['Area'] * one_pix)
                                           
        cell_df.drop(columns='Area', inplace=True)
        
    adata = ad.AnnData(cell_df)
    
    adata.layers['raw'] = np.copy(adata.X)
    
    for key in obs_anno.keys():
        adata.obs[key] = obs_anno[key]
        
    return adata

# %%

cell_dfs = []

for filepath in filepaths:
    cell_dfs.append(load_cell_df(filepath, usecols, index_col='Cell ID'))
    
cell_df = concat_cell_dfs(cell_dfs)

adata = build_cell_anndata(cell_df, obs_anno, norm_by_area=False)
    
# %%

sc.pp.log1p(adata)
sc.pp.scale(adata)
sc.tl.pca(adata)

sc.pp.neighbors(adata)
sc.tl.umap(adata, min_dist=0.5)
sc.tl.leiden(adata, flavor="igraph", n_iterations=-1)

# %%

sc.pl.umap(adata, color='leiden')

# %%

from sklearn.preprocessing import minmax_scale, normalize

raw_exp = np.copy(adata.layers['raw'])
# raw_exp = normalize(raw_exp, norm='l1')
# raw_exp = normalize(raw_exp, norm='l1')
raw_exp = raw_exp[:,3:]
area = raw_exp[:,-1]
area = np.expand_dims(area, axis=1)
raw_exp = np.delete(raw_exp, -1, axis=1)
# raw_exp[:,:3] = raw_exp[:,:3] * area
raw_exp = raw_exp / area
# raw_exp = minmax_scale(raw_exp)

R = np.corrcoef(np.transpose(raw_exp))

fig, ax = plt.subplots()
im = ax.imshow(R, vmin=-0.5, vmax=0.5)
fig.colorbar(im, ax=ax)

genes = adata.var_names[:-1]
genes = genes[3:]
ax.set_xticks(range(len(genes)), labels=genes,
              rotation=45, ha="right", rotation_mode="anchor")
ax.set_yticks(range(len(genes)), labels=genes)

# for i in range(len(genes)):
#     for j in range(len(genes)):
#         text = ax.text(j, i, np.round(R[i, j], decimals=1),
#                        ha="center", va="center", color="w")
        
plt.title('HCR data correlation matrix (normalized by mask area)')

# %%

import pingouin as pg

raw_exp = np.copy(adata.layers['raw'])
# raw_exp = raw_exp[:,3:-1]
# genes = adata.var_names[3:-1]
raw_exp = raw_exp[:,:-1]
genes = adata.var_names[:-1]
df = pd.DataFrame(raw_exp, columns=genes)

R = pg.pcorr(df)

fig, ax = plt.subplots()
im = ax.imshow(R, vmin=-0.5, vmax=0.5)
fig.colorbar(im, ax=ax)

ax.set_xticks(range(len(genes)), labels=genes,
              rotation=45, ha="right", rotation_mode="anchor")
ax.set_yticks(range(len(genes)), labels=genes)

plt.title('HCR data correlation matrix (partial correlation, nonbinary genes)')

#%%

import pandas as pd
import numpy as np
import pingouin as pg
from sklearn.linear_model import LinearRegression

def partial_pointbiserial_matrix(df):
    """
    Computes partial point-biserial correlations between all binary and
    continuous variables in df, controlling for all other continuous variables
    (matching the covariate structure of pg.pcorr() on the continuous block).

    Binary variables are identified by '_b' in their column name.
    All other columns are treated as continuous.

    For each continuous variable, residualizes on all OTHER continuous
    variables, matching the pairwise partial correlation convention.

    Parameters
    ----------
    df : pd.DataFrame
        No NaNs. Column names containing '_b' are treated as binary and must
        have exactly 2 unique values. All other columns are treated as
        continuous.

    Returns
    -------
    r_matrix : pd.DataFrame
        Shape (n_binary, n_continuous), correlation coefficients.

    p_matrix : pd.DataFrame
        Shape (n_binary, n_continuous), two-tailed p-values.

    results_long : pd.DataFrame
        Tidy table with columns: binary_var, continuous_var, n, r, CI95%, p-val.

    Example
    -------
        df = pd.DataFrame({
            'treatment_b':   [0, 1, 0, 1, 0, 1, 0, 1],
            'female_b':      [1, 1, 0, 0, 1, 0, 1, 0],
            'smoker_b':      [0, 0, 1, 0, 1, 1, 0, 0],
            'reaction_time': [520, 480, 610, 455, 530, 490, 580, 460],
            'memory_score':  [74, 88, 65, 91, 70, 85, 68, 93],
            'anxiety_score': [12, 8, 15, 6, 11, 9, 14, 7],
        })

        r_matrix, p_matrix, results_long = partial_pointbiserial_matrix(df)
    """

    binary_cols     = [c for c in df.columns if '_b' in c]
    continuous_cols = [c for c in df.columns if '_b' not in c]

    # --- input validation ---
    if df.isnull().any().any():
        raise ValueError("NaNs detected. Drop or impute before calling this function.")

    for col in binary_cols:
        if df[col].nunique() != 2:
            raise ValueError(f"'{col}' must have exactly 2 unique values.")

    # --- residualize each continuous variable on all OTHER continuous variables ---
    residuals = {}
    for cont_col in continuous_cols:
        other_cols = [c for c in continuous_cols if c != cont_col]
        X = df[other_cols].values
        y = df[cont_col].values
        residuals[cont_col] = y - LinearRegression().fit(X, y).predict(X)

    # --- point-biserial for every binary-continuous pair ---
    rows = []
    for bin_col in binary_cols:
        for cont_col in continuous_cols:
            temp = pd.DataFrame({
                bin_col:     df[bin_col].values,
                'residuals': residuals[cont_col]
            })
            result = pg.corr(temp[bin_col], temp['residuals'], method='pearson')
            rows.append({
                'binary_var':     bin_col,
                'continuous_var': cont_col,
                'n':              result['n'].iloc[0],
                'r':              result['r'].iloc[0],
                'CI95':          result['CI95'].iloc[0],
                'p_val':          result['p_val'].iloc[0],
            })

    results_long = pd.DataFrame(rows)

    r_matrix = results_long.pivot(
        index='binary_var', columns='continuous_var', values='r'
    ).loc[binary_cols, continuous_cols]

    p_matrix = results_long.pivot(
        index='binary_var', columns='continuous_var', values='p_val'
    ).loc[binary_cols, continuous_cols]

    return r_matrix, p_matrix, results_long

r_matrix, p_matrix, results_long = partial_pointbiserial_matrix(df)

binary_cols     = [c for c in df.columns if '_b' in c]
continuous_cols = [c for c in df.columns if '_b' not in c]
R.loc[binary_cols, continuous_cols] = r_matrix.values
R.loc[continuous_cols, binary_cols] = r_matrix.T.values

fig, ax = plt.subplots()
im = ax.imshow(R, vmin=-0.5, vmax=0.5)
cbar = fig.colorbar(im, ax=ax)
cbar.set_label('Partial Correlation Coefficient', rotation=270, fontsize=12)

ax.set_xticks(range(len(genes)), labels=genes,
              rotation=45, ha="right", rotation_mode="anchor")
ax.set_yticks(range(len(genes)), labels=genes)
        
plt.title('HCR data correlation matrix (partial correlation)')


# %%

import os
from ABC_toolbox import cell_funcs, gene_funcs

n_real_cells = len(cell_df)

local_data_dir = r'C:\Users\jpv88\Documents\GenePicker9001_data'
abc_genes = np.load(r'C:\Users\jpv88\Documents\GitHub\GenePicker9001\ABC_toolbox\util_files\gene_scRNAseq.npy', 
                    allow_pickle=True)
abc_genes = [gene.lower() for gene in abc_genes]
genes = ['slc17a6', 'slc32a1', 'slc5a7', 'phox2b', 'ralyl', 'tenm2', 'ebf3', 'pcp4', 'tshz2',
         'alcam', 'celf2', 'meis2', 'rph3a', 'robo1', 'syt1', 'zfhx3']
abc_genes_idx = [np.where(gene == np.array(abc_genes))[0][0] for gene in genes]

meta = pd.read_csv(os.path.join(local_data_dir, "antIRN-PARN-scRNAseq-meta.csv"), low_memory=False)
exp = np.load(os.path.join(local_data_dir, "antIRN-PARN-scRNAseq-raw.npy"))
freqs = pd.read_pickle(os.path.join(local_data_dir, "antIRN-PARN-MERFISH-freqs.pkl"))
exp = exp[:,abc_genes_idx]

exp_super, meta_super = cell_funcs.boot_super(exp, meta, k=1)

# bootstrap distribution from scRNAseq that matches MERFISH frequencies
exp_boot, meta_boot = cell_funcs.bootstrap_scRNAseq(meta_super, exp_super, freqs, 
                                                    n=n_real_cells)

# %%

df = pd.DataFrame(exp_boot, columns=genes)
R = pg.pcorr(df)
# R = np.corrcoef(np.transpose(exp_boot))

fig, ax = plt.subplots()
im = ax.imshow(R, vmin=-0.5, vmax=0.5)
fig.colorbar(im, ax=ax)

# genes = adata.var_names
# genes = genes[0:-1]
ax.set_xticks(range(len(genes)), labels=genes,
              rotation=45, ha="right", rotation_mode="anchor")
ax.set_yticks(range(len(genes)), labels=genes)

# for i in range(len(genes)):
#     for j in range(len(genes)):
#         text = ax.text(j, i, np.round(R[i, j], decimals=1),
#                        ha="center", va="center", color="w")
        
plt.title('scRNA-seq (ABC atlas) correlation matrix')

#%% quantifying gene expression im 

import tifffile

from analysis_utils import matrix_utils, filter_data, img_utils
from tqdm import tqdm

data_dir = 'D:\\2026-01-16_MC_SC_17\\post\\core_output\\s03L'

uniq_id, post_dir, exp_name = matrix_utils.derive_params(data_dir)
data_f, reg_f, bin_f, mask_f, mask_im_f = matrix_utils.build_paths(data_dir, 
                                                                   uniq_id, 
                                                                   post_dir, 
                                                                   exp_name,
                                                                   mode='curated')

cell_dfs = []

for filepath in filepaths:
    cell_dfs.append(load_cell_df(filepath, None, index_col='Cell ID'))
    
cell_df = concat_cell_dfs(cell_dfs)

adata = build_cell_anndata(cell_df, obs_anno, norm_by_area=False)

masks_im = tifffile.imread(mask_im_f)
masks_im = np.transpose(masks_im,[1,2,0])

pp_masks_im = filter_data.filt_masks_im_pp(masks_im, cell_df)

def color_masks_by_gene(masks_im, cell_df, gene):
    
    gene_exp = cell_df[gene].values
    gene_exp = gene_exp.astype(np.uint16)
    masks = cell_df['Mask ID'].values
    
    lut = np.zeros(masks.max() + 1, dtype=np.uint16)
    lut[masks] = gene_exp
    
    gene_exp_im = lut[masks_im]
    
    return gene_exp_im

gene_exp_im = color_masks_by_gene(pp_masks_im, cell_df, 'phox2b')
tifffile.imwrite('phox2b_exp.tif', gene_exp_im.transpose(2, 0, 1))
# %%

for gene in cell_df.columns:

    sc.pl.umap(adata, color=gene)
    
# %%

sc.pl.umap(adata, color='egfp')

# %%

sc.pl.pca(adata, color='slc17a6')

