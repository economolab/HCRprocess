# -*- coding: utf-8 -*-
"""

"""

import anndata as ad
import numpy as np
import pandas as pd

from ABC_toolbox import ABC_utils, cell_funcs

import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.rcParams['image.composite_image'] = False
plt.rcParams['svg.fonttype'] = 'none'

# %%

filepaths = ['MC_SC_17_s03L.csv', 'MC_SC_17_s01R.csv']

usecols = ['Cell ID', 'slc17a6', 'slc32a1', 'slc5a7', 'phox2b', 
           'ralyl', 'tenm2', 'ebf3', 'pcp4', 'tshz2', 'alcam', 'celf2', 
           'meis2', 'rph3a', 'robo1', 'syt1', 'zfhx3']

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

HCR_exp = cell_df.to_numpy()
HCR_genes = cell_df.columns
HCR_genes = [gene.capitalize() for gene in HCR_genes]

cht_group = ['242 PGRNd Dmbx1 Glut', '252 DMX VII Tbx20 Chol', '261 HB Calcb Chol']
vgat_group = ['264 PRNc Otp Gly-Gaba','266 PRNc Prox1 Brs3 Gly-Gaba','278 NLL Gata3 Gly-Gaba',
              '282 POR Spp1 Gly-Gaba','283 PRP Otp Gly-Gaba','285 MY Lhx1 Gly-Gaba',
              '287 MV-SPIV-PRP Dmbx1 Gly-Gaba','290 MY Prox1 Lmo7 Gly-Gaba','299 MARN-PPY Ngfr Gly-Gaba',
              '300 PARN-MDRNd-NTS Gbx2 Gly-Gaba','301 MV Nr4a2 Gly-Gaba','302 MV Xdh Gly-Gaba',
              '303 IRN Dmbx1 Pax2 Gly-Gaba','304 NTS-PARN Neurod2 Gly-Gaba','306 SPVI-SPVC Sall3 Lhx1 Gly-Gaba',
              '308 DCO Il22 Gly-Gaba']
vglut_group = ['224 PCG-PRNr Vsx2 Nkx6-1 Glut','225 PRNc-NI-SG-RPO Vsx2 Nr4a2 Glut',
               '226 PRNc-PARN Tlx1 Glut','228 PSV Pvalb Lhx2 Glut','229 PB-NTS Phox2b Ebf3 Lmx1b Glut',
               '230 PRNr Otp Nfib Glut','233 NLL-SOC Spp1 Glut','235 PG-TRN-LRN Fat2 Glut',
               '236 IRN Vip Glut','237 PRP-NI-PRNc-GRN Otp Glut','238 NTS Phox2b Glut','243 PGRN-PARN-MDRN Hoxb5 Glut',
               '245 SPVI-SPVC Tlx3 Ebf3 Glut','246 CU-ECU-SPVI Foxb1 Glut','247 MV-SPIV Phox2b Ebf3 Lbx1 Glut',
               '254 VCO Mafa Meis2 Glut','255 SPVO Mafa Meis2 Glut']

cell_df['nt_group'] = np.repeat('None', len(cell_df))
cht_bool = (cell_df['slc5a7'] == 1)
cell_df.loc[cht_bool, 'nt_group'] = 'cht'
vgat_bool = (cell_df['slc32a1'] == 1) & (cell_df['nt_group'] == 'None')
cell_df.loc[vgat_bool, 'nt_group'] = 'vgat'
vglut_bool = (cell_df['slc17a6'] == 1) & (cell_df['nt_group'] == 'None')
cell_df.loc[vglut_bool, 'nt_group'] = 'vglut'




#%%

import os

n_real_cells = len(cell_df)

local_data_dir = r'C:\Users\jpv88\Documents\GenePicker9001_data'
meta = pd.read_csv(os.path.join(local_data_dir, "antIRN-PARN-scRNAseq-meta.csv"), low_memory=False)
exp = np.load(os.path.join(local_data_dir, "antIRN-PARN-scRNAseq-norm.npy"))
freqs = pd.read_pickle(os.path.join(local_data_dir, "antIRN-PARN-MERFISH-freqs.pkl"))

exp_super, meta_super = cell_funcs.boot_super(exp, meta, k=3)
boot_mat, boot_meta = cell_funcs.bootstrap_scRNAseq(meta_super, exp_super, freqs, n=n_real_cells)

ABC_genes = list(ABC_utils.load_gene('scRNAseq'))
ABC_gene_idx = [ABC_genes.index(gene) for gene in HCR_genes]

boot_mat = boot_mat[:,ABC_gene_idx]

# %%

from anndata import AnnData

adata = AnnData(np.vstack((HCR_exp, boot_mat)))

# assign batch
n_HCR = HCR_exp.shape[0]
n_ABC = boot_mat.shape[0]

adata.obs['batch'] = np.concatenate((n_HCR*['HCR'], n_ABC*['ABC']))
adata.var_names = HCR_genes

# %%

import scvi
import scanpy.external as sce

n_layers = 3
n_latent = 16
n_neighbors = 40
latent_key = "X_pca_harmony"
    
scvi.model.SCVI.setup_anndata(adata, batch_key="batch")

model = scvi.model.SCVI(
    adata,
    n_layers=n_layers,
    n_latent=n_latent
)

model.train(max_epochs=100,
            early_stopping=True,
            early_stopping_monitor="reconstruction_loss_validation",
            early_stopping_patience=2)

adata.obsm["X_scVI"] = model.get_latent_representation(adata)

sce.pp.harmony_integrate(adata, 
                         "batch",
                         basis="X_scVI",
                         max_iter_harmony=100,
                         max_iter_kmeans=20)

#%%

level = 'subclass'

HCR_bool = (adata.obs['batch'].values == 'HCR')
ABC_bool = (adata.obs['batch'].values == 'ABC')

from sklearn.neighbors import KNeighborsClassifier
neigh = KNeighborsClassifier(n_neighbors=5)

X = adata.obsm['X_pca_harmony']
X = X[ABC_bool,:]
y = boot_meta[level].values
neigh.fit(X, y)

X = adata.obsm['X_pca_harmony']
X = X[HCR_bool,:]
y = neigh.predict(X)

cell_df[level] = y

#%%

def calc_frac_per_type(cell_labels):
    
    # find every unique cell type of the prescribed level in this metadata dataframe
    uniq_clu = np.unique(cell_labels)
    
    # find how many cells of each type there are
    num_per_clu = []
    for clu in uniq_clu:
        num_per_clu.append(sum(cell_labels == clu))
    
    # calculate the fraction per cell type
    num_cells = len(cell_labels)
    frac_per_clu = [x/num_cells for x in num_per_clu]
    
    # convert to dictionary
    d = {uniq_clu[i]: frac_per_clu[i] for i in range(len(uniq_clu))}
    
    return d

if level == 'cluster':
    subclass_freqs = freqs
else:
    subclass_freqs = cell_funcs.recalc_freqs(freqs, new_level=level)
    
pred_freqs = calc_frac_per_type(y)

#%%

# Sample data
real_subclasses = np.unique(list(subclass_freqs.keys()))
categories = real_subclasses
group_a = list(subclass_freqs.values())
group_b = []

for real_subclass in real_subclasses:
    if real_subclass in pred_freqs:
        group_b.append(pred_freqs[real_subclass])
    else:
        group_b.append(0)
        
group_a = np.array(group_a) * 100
group_b = np.array(group_b) * 100
        
x = np.arange(len(categories))
width = 0.35  # width of each bar
 
fig, ax = plt.subplots(figsize=(8, 5))
 
bars_a = ax.bar(x - width/2, group_a, width, label='ABC Atlas Cell Type Proportions', color='steelblue')
bars_b = ax.bar(x + width/2, group_b, width, label='HCR Data Cell Type Proportions', color='coral')

# Labels and formatting
ax.set_xlabel('ABC Atlas ant-IRN-PARN Subclasses')
ax.set_ylabel('Cell Type Proportion (%)')
ax.set_title('ABC Atlas vs. HCR Data')
ax.legend()

# %%

usecols = ['Cell ID', 'Mask ID', 'slc17a6', 'slc32a1', 'slc5a7', 'phox2b', 
           'ralyl', 'tenm2', 'ebf3', 'pcp4', 'tshz2', 'alcam', 'celf2', 
           'meis2', 'rph3a', 'robo1', 'syt1', 'zfhx3']

cell_dfs = []

for filepath in filepaths:
    cell_dfs.append(load_cell_df(filepath, usecols, index_col='Cell ID'))
    
cell_df = concat_cell_dfs(cell_dfs)
cell_df[level] = y

import tifffile
from ABC_toolbox import ABC_plot
from tqdm import tqdm

colors_dict = ABC_plot.fetch_colors_dict('subclass')

for key, value in colors_dict.items():
    new_value = np.round(np.array(value) * 255)
    new_value = new_value.astype(np.uint8)
    colors_dict[key] = new_value
    
masks_im = tifffile.imread('MC_SC_17_s03L_pp_masks.tif')

z, h, w = np.shape(masks_im)
rgb_masks = np.zeros((z, h, w, 4), dtype=np.uint8)

in_im_bool = np.array(['s03L' in cell_id for cell_id in cell_df.index.values])
in_im_cell_df = cell_df[in_im_bool]

for i in tqdm(range(len(in_im_cell_df))):
    mask_id = in_im_cell_df['Mask ID'].iloc[i]
    im_mask = (masks_im == mask_id)
    c = colors_dict[in_im_cell_df['subclass'].iloc[i]]
    c = np.append(c, 255)
    c = c.astype(np.uint8)
    c = c.reshape(1, -1)
    rgb_masks[im_mask,:] = c

tifffile.imwrite(
    's03L_subclasses.tiff',
    rgb_masks,
    photometric='rgb',
    extrasamples=['unassalpha'],
    imagej=True
)
#%%

norm_exp = model.get_normalized_expression(adata,transform_batch='HCR')

#%%



def plot_hist_gene(adata, gene):
    
    norm_cell_df = model.get_normalized_expression(adata,transform_batch='HCR')
    gene_exp = norm_cell_df[gene].values
    
    HCR_bool = (adata.obs['batch'].values == 'HCR')
    ABC_bool = (adata.obs['batch'].values == 'ABC')
    
    plt.hist(gene_exp[HCR_bool], bins=20, alpha=0.5, label='HCR')
    plt.hist(gene_exp[ABC_bool], bins=20, alpha=0.5, label='ABC')
    plt.legend()
    
plot_hist_gene(adata,'')

# %%

import os

n_real_cells = len(cell_df)

local_data_dir = r'C:\Users\jpv88\Documents\GenePicker9001_data'
meta = pd.read_csv(os.path.join(local_data_dir, "antIRN-PARN-scRNAseq-meta.csv"), low_memory=False)
exp = np.load(os.path.join(local_data_dir, "antIRN-PARN-scRNAseq-norm.npy"))
freqs = pd.read_pickle(os.path.join(local_data_dir, "antIRN-PARN-MERFISH-freqs.pkl"))

exp_super, meta_super = cell_funcs.boot_super(exp, meta, k=3)
boot_mat, boot_meta = cell_funcs.bootstrap_scRNAseq(meta_super, exp_super, freqs, n=n_real_cells)

ABC_plot.plot_gene('Slc17a6', meta_super, exp_super, level='subclass')
ABC_plot.plot_gene('Slc5a7', meta_super, exp_super, level='subclass')
ABC_plot.plot_gene('Slc32a1', meta_super, exp_super, level='subclass')
ABC_plot.plot_gene('Phox2b', meta_super, exp_super, level='subclass')
    