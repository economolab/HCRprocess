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
           'meis2', 'rph3a', 'robo1', 'syt1', 'zfhx3']

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

#%%

# 1. Scale only the continuous features, leave binary as-is (0/1)
continuous_cols = ['phox2b', 'ralyl', 'tenm2', 'ebf3', 'pcp4', 'tshz2', 
                   'alcam', 'celf2', 'meis2', 'rph3a', 'robo1', 'syt1', 'zfhx3']
binary_cols = ['slc17a6_b', 'slc32a1_b', 'slc5a7_b']

adata[:, continuous_cols].X = sc.pp.scale(adata[:, continuous_cols], copy=True).X

# 2. Downweight binary columns
for col in binary_cols:
    idx = adata.var_names.get_loc(col)
    adata.X[:, idx] *= 2  # tune this

# 3. Skip PCA, run neighbors directly on the 16-dim weighted feature space
sc.pp.neighbors(adata, n_neighbors=10, use_rep='X')
sc.tl.umap(adata, min_dist=0.5)
sc.tl.leiden(adata, flavor="igraph", n_iterations=-1)
    
# %%

# sc.pp.log1p(adata)
sc.pp.scale(adata)

binary_col_idx = [0, 1, 2]
X = adata.X.copy()
weight = 0.2
for col in binary_col_idx:
    X[:, col] *= weight
    
sc.tl.pca(adata)

sc.pp.neighbors(adata)
sc.tl.umap(adata, min_dist=0.5)
# sc.tl.leiden(adata, flavor="igraph", n_iterations=-1)

# %%
# sc.pl.umap(adata)
sc.pl.umap(adata, color='leiden')

#%%
sc.tl.pca(adata, n_comps=2)
sc.pl.pca(adata, color='leiden')

# %%

import pingouin as pg

raw_exp = np.copy(adata.layers['raw'])
genes = adata.var_names
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
    binary_cols     = [c for c in df.columns if '_b' in c]
    continuous_cols = [c for c in df.columns if '_b' not in c]

    if df.isnull().any().any():
        raise ValueError("NaNs detected. Drop or impute before calling this function.")

    for col in binary_cols:
        if df[col].nunique() != 2:
            raise ValueError(f"'{col}' must have exactly 2 unique values.")

    # residualize each continuous variable on all OTHER continuous variables
    cont_residuals = {}
    for cont_col in continuous_cols:
        other_cols = [c for c in continuous_cols if c != cont_col]
        X = df[other_cols].values
        y = df[cont_col].values
        cont_residuals[cont_col] = y - LinearRegression().fit(X, y).predict(X)

    # residualize each binary variable on the SAME covariate set used for
    # each continuous variable (i.e. all continuous vars EXCEPT cont_col) —
    # this must be recomputed per (bin_col, cont_col) pair, not once per bin_col
    bin_residuals = {}
    for cont_col in continuous_cols:
        other_cols = [c for c in continuous_cols if c != cont_col]
        X = df[other_cols].values
        for bin_col in binary_cols:
            y = df[bin_col].values.astype(float)
            bin_residuals[(bin_col, cont_col)] = y - LinearRegression().fit(X, y).predict(X)

    rows = []
    n_covar = len(continuous_cols) - 1
    n = len(df)
    dof = n - 2 - n_covar

    from scipy import stats
    for bin_col in binary_cols:
        for cont_col in continuous_cols:
            x = bin_residuals[(bin_col, cont_col)]
            y = cont_residuals[cont_col]

            r = np.corrcoef(x, y)[0, 1]
            t_stat = r * np.sqrt(dof / (1 - r**2))
            p_corrected = 2 * stats.t.sf(np.abs(t_stat), dof)

            rows.append({
                'binary_var':     bin_col,
                'continuous_var': cont_col,
                'n':              n,
                'r':              r,
                'p_val':          p_corrected,
            })

    results_long = pd.DataFrame(rows)
    r_matrix = results_long.pivot(index='binary_var', columns='continuous_var', values='r').loc[binary_cols, continuous_cols]
    p_matrix = results_long.pivot(index='binary_var', columns='continuous_var', values='p_val').loc[binary_cols, continuous_cols]

    return r_matrix, p_matrix, results_long

genes = adata.var_names
raw_exp = np.copy(adata.layers['raw'])
df = pd.DataFrame(raw_exp, columns=genes)
r_matrix, p_matrix, results_long = partial_pointbiserial_matrix(df)

binary_cols     = [c for c in df.columns if '_b' in c]
continuous_cols = [c for c in df.columns if '_b' not in c]
R.loc[binary_cols, continuous_cols] = r_matrix.values
R.loc[continuous_cols, binary_cols] = r_matrix.T.values

fig, ax = plt.subplots()

cmap = plt.cm.viridis.copy()
cmap.set_bad(color='gainsboro')
R_masked = np.ma.masked_where(np.eye(R.shape[0], dtype=bool), R.values)
corr_HCR = R_masked.copy()

im = ax.imshow(R_masked, vmin=-0.5, vmax=0.5)
cbar = fig.colorbar(im, ax=ax)
cbar.set_label('Pearson Partial Correlation Coefficient', rotation=270, fontsize=12, labelpad=20)

ax.set_xticks(range(len(genes)), labels=genes,
              rotation=45, ha="right", rotation_mode="anchor")
ax.set_yticks(range(len(genes)), labels=genes)
        
plt.title('HCR data correlation matrix (partial correlation)')




# %%

import os
from ABC_toolbox import cell_funcs

n_real_cells = len(cell_df)

local_data_dir = r'C:\Users\jpv88\Documents\GenePicker9001_data'
abc_genes = np.load(r'C:\Users\jpv88\Documents\GitHub\GenePicker9001\ABC_toolbox\util_files\gene_scRNAseq.npy', 
                    allow_pickle=True)
abc_genes = [gene.lower() for gene in abc_genes]
genes = ['slc17a6', 'slc32a1', 'slc5a7', 'phox2b', 'ralyl', 'tenm2', 'ebf3', 'pcp4', 'tshz2',
         'alcam', 'celf2', 'meis2', 'rph3a', 'robo1', 'syt1', 'zfhx3']
abc_genes_idx = [np.where(gene == np.array(abc_genes))[0][0] for gene in genes]

meta = pd.read_csv(os.path.join(local_data_dir, "antIRN-PARN-scRNAseq-meta.csv"), low_memory=False)
exp = np.load(os.path.join(local_data_dir, "antIRN-PARN-scRNAseq-norm.npy"))
freqs = pd.read_pickle(os.path.join(local_data_dir, "antIRN-PARN-MERFISH-freqs.pkl"))
exp = exp[:,abc_genes_idx]

exp_super, meta_super = cell_funcs.boot_super(exp, meta, k=1)

# bootstrap distribution from scRNAseq that matches MERFISH frequencies
exp_boot, meta_boot = cell_funcs.bootstrap_scRNAseq(meta_super, exp_super, freqs, 
                                                    n=n_real_cells)

df_ABC = pd.DataFrame(exp_boot, columns=genes)
R = pg.pcorr(df_ABC)
np.fill_diagonal(R.values, 0)

cmap = plt.cm.viridis.copy()
cmap.set_bad(color='gainsboro')

R_masked = np.ma.masked_where(np.eye(R.shape[0], dtype=bool), R.values)
corr_ABC = R_masked.copy()

fig, ax = plt.subplots()
im = ax.imshow(R_masked, cmap=cmap, vmin=-0.5, vmax=0.5)
cbar = fig.colorbar(im, ax=ax)
cbar.set_label('Pearson Partial Correlation Coefficient', rotation=270, fontsize=12, labelpad=20)
ax.set_xticks(range(len(genes)), labels=genes, rotation=45, ha="right", rotation_mode="anchor")
ax.set_yticks(range(len(genes)), labels=genes)
plt.title("scRNA-seq (ABC atlas) correlation matrix")


# %%

for gene in cell_df.columns:

    sc.pl.umap(adata, color=gene)
    
# %%

sc.pl.umap(adata, color='egfp')

# %%

sc.pl.pca(adata, color='slc17a6')

#%%

spearman_corr = cell_df.corr(method='spearman')
corr_masked = spearman_corr.copy()
np.fill_diagonal(corr_masked.values, np.nan)

fig, ax = plt.subplots(figsize=(8, 6))
im = ax.imshow(corr_masked, cmap='coolwarm', vmin=-1, vmax=1)
ax.set_xticks(range(len(corr_masked.columns)))
ax.set_yticks(range(len(corr_masked.columns)))
ax.set_xticklabels(corr_masked.columns, rotation=45, ha='right')
ax.set_yticklabels(corr_masked.columns)
fig.colorbar(im)
plt.tight_layout()
plt.show()

#%%

import matplotlib.pyplot as plt

# ---- global style, set once outside the loop ----
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'Helvetica', 'DejaVu Sans'],
    'axes.titlesize': 13,
    'axes.titleweight': 'bold',
    'axes.labelsize': 11,
    'xtick.labelsize': 9,
    'ytick.labelsize': 9,
    'axes.edgecolor': '#333333',
    'figure.dpi': 120,
})

colors = ['#4C72B0', '#DD8452']

for gene in genes:

    fig, ax = plt.subplots(figsize=(6, 4))

    data_HCR = np.asarray(adata[:, gene].X).flatten()
    data_ABC = df_ABC[gene].values

    bins = np.linspace(min(data_HCR.min(), data_ABC.min()),
                        max(data_HCR.max(), data_ABC.max()), 30)

    ax.hist(data_HCR, bins=bins, alpha=0.55, label=f'HCR (n={len(data_HCR)})',
            density=True, color=colors[0], linewidth=0.6)
    ax.hist(data_ABC, bins=bins, alpha=0.55, label=f'ABC (n={len(data_ABC)})',
            density=True, color=colors[1], linewidth=0.6)

    # mean markers as small dashed vertical lines
    for data, color in zip([data_HCR, data_ABC], colors):
        ax.axvline(data.mean(), color=color, linestyle='--', linewidth=1, alpha=0.8)

    ax.set_title(gene, pad=10)
    ax.set_xlabel('Counts')
    ax.set_ylabel('Probability Density')

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(0.8)
    ax.spines['bottom'].set_linewidth(0.8)

    ax.grid(axis='y', alpha=0.25, linewidth=0.6)
    ax.set_axisbelow(True)

    ax.legend(frameon=False, loc='upper right')

    fig.tight_layout()
    
#%%

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

def plot_split_corr(corr_a, corr_b, labels=None, title_a="A", title_b="B",
                     cmap="viridis", vmin=-0.5, vmax=0.5, figsize=(8, 8)):
    """
    Overlay two correlation matrices: upper triangle = corr_a, lower triangle = corr_b.
    corr_a, corr_b: square pandas DataFrames or numpy arrays, same shape/order.
    """
    a = corr_a.values if isinstance(corr_a, pd.DataFrame) else np.asarray(corr_a)
    b = corr_b.values if isinstance(corr_b, pd.DataFrame) else np.asarray(corr_b)
    n = a.shape[0]
    if labels is None:
        labels = corr_a.columns if isinstance(corr_a, pd.DataFrame) else [str(i) for i in range(n)]

    combined = np.full((n, n), np.nan)
    iu = np.triu_indices(n, k=1)   # upper triangle (excl diagonal)
    il = np.tril_indices(n, k=-1)  # lower triangle (excl diagonal)
    combined[iu] = a[iu]
    combined[il] = b[il]

    masked = np.ma.masked_invalid(combined)

    cmap_obj = plt.get_cmap(cmap).copy()
    cmap_obj.set_bad(color="lightgray")  # diagonal (masked/NaN) cells

    fig, ax = plt.subplots(figsize=figsize)
    norm = Normalize(vmin=vmin, vmax=vmax)
    im = ax.imshow(masked, cmap=cmap_obj, norm=norm, interpolation="nearest")

    ax.set_xticks(range(n)); ax.set_xticklabels(labels, rotation=45, ha="right")
    ax.set_yticks(range(n)); ax.set_yticklabels(labels)
    ax.set_title(f"Upper: {title_a}   |   Lower: {title_b}")

    fig.colorbar(im, ax=ax, shrink=0.8, label="correlation")
    plt.tight_layout()
    return fig, ax

diff = np.abs(corr_HCR - corr_ABC)
plt.imshow(diff, cmap="viridis", vmin=0, vmax=1)
plt.colorbar(label="corr_a - corr_b")

#%%

from scipy.stats import rankdata

# --- continuous-continuous block: full matrix via pg.pcorr on ranks ---
binary_cols     = [c for c in df.columns if '_b' in c]
continuous_cols = [c for c in df.columns if '_b' not in c]
df_ranked_full = df.copy()
df_ranked_full[continuous_cols] = df[continuous_cols].apply(rankdata)
R = pg.pcorr(df_ranked_full)   # NOTE: this now includes binary cols un-ranked,
                                # matching your original mixed-matrix approach

# --- binary-continuous block ---
def partial_pointbiserial_matrix_spearman(df):
    binary_cols     = [c for c in df.columns if '_b' in c]
    continuous_cols = [c for c in df.columns if '_b' not in c]

    if df.isnull().any().any():
        raise ValueError("NaNs detected. Drop or impute before calling this function.")

    for col in binary_cols:
        if df[col].nunique() != 2:
            raise ValueError(f"'{col}' must have exactly 2 unique values.")

    # rank-transform continuous variables once
    ranks = {c: rankdata(df[c].values) for c in continuous_cols}

    # residualize ranked cont_col on ranks of all OTHER continuous variables
    cont_residuals = {}
    for cont_col in continuous_cols:
        other_cols = [c for c in continuous_cols if c != cont_col]
        X = np.column_stack([ranks[c] for c in other_cols])
        y = ranks[cont_col]
        cont_residuals[cont_col] = y - LinearRegression().fit(X, y).predict(X)

    # residualize binary vars on the SAME ranked covariate set, per cont_col
    bin_residuals = {}
    for cont_col in continuous_cols:
        other_cols = [c for c in continuous_cols if c != cont_col]
        X = np.column_stack([ranks[c] for c in other_cols])
        for bin_col in binary_cols:
            y = df[bin_col].values.astype(float)
            bin_residuals[(bin_col, cont_col)] = y - LinearRegression().fit(X, y).predict(X)

    rows = []
    n_covar = len(continuous_cols) - 1
    n = len(df)
    dof = n - 2 - n_covar

    from scipy import stats
    for bin_col in binary_cols:
        for cont_col in continuous_cols:
            x = bin_residuals[(bin_col, cont_col)]
            y = cont_residuals[cont_col]

            r = np.corrcoef(x, y)[0, 1]
            t_stat = r * np.sqrt(dof / (1 - r**2))
            p_corrected = 2 * stats.t.sf(np.abs(t_stat), dof)

            rows.append({
                'binary_var':     bin_col,
                'continuous_var': cont_col,
                'n':              n,
                'r':              r,
                'p_val':          p_corrected,
            })

    results_long = pd.DataFrame(rows)
    r_matrix = results_long.pivot(index='binary_var', columns='continuous_var', values='r').loc[binary_cols, continuous_cols]
    p_matrix = results_long.pivot(index='binary_var', columns='continuous_var', values='p_val').loc[binary_cols, continuous_cols]

    return r_matrix, p_matrix, results_long

r_matrix, p_matrix, results_long = partial_pointbiserial_matrix_spearman(df)
R.loc[binary_cols, continuous_cols] = r_matrix.values
R.loc[continuous_cols, binary_cols] = r_matrix.T.values

fig, ax = plt.subplots()

cmap = plt.cm.viridis.copy()
cmap.set_bad(color='gainsboro')
R_masked = np.ma.masked_where(np.eye(R.shape[0], dtype=bool), R.values)

im = ax.imshow(R_masked, vmin=-0.5, vmax=0.5)
cbar = fig.colorbar(im, ax=ax)
cbar.set_label('Spearman Partial Correlation Coefficient', rotation=270, fontsize=12)

ax.set_xticks(range(len(genes)), labels=genes,
              rotation=45, ha="right", rotation_mode="anchor")
ax.set_yticks(range(len(genes)), labels=genes)
        
plt.title('HCR data correlation matrix (partial correlation)')

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

#%%

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


