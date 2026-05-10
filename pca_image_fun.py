# -*- coding: utf-8 -*-
"""

"""

#%%

import tifffile

import numpy as np

from sklearn.preprocessing import StandardScaler, MinMaxScaler
from sklearn.decomposition import PCA
from tqdm import tqdm

from analysis_utils import load_data, filter_data, img_utils

#%%

exp_dirs = [r'D:\2026-01-16_MC_SC_17']
masks_im_f = r'D:\2026-01-16_MC_SC_17\post\masks\s03L\s03L_masks_qc_final.tif'
neur_im_f = r'D:\2026-01-16_MC_SC_17\post\core_output\s03L\r1_neur445.tif'

n_components = 6
min_pc_percentile = 10
max_pc_percentile = 90

#%%

cell_df_orig = load_data.load_experiments(exp_dirs)
cell_df_orig = filter_data.filt_cells(cell_df_orig, 's03L')

cell_df_gene = load_data.load_experiments(exp_dirs)
cell_df_gene = filter_data.filt_cells(cell_df_gene, 's03L')
cell_df_gene = filter_data.filt_features(cell_df_gene, remove_features=['dtom', 'egfp'], genes=True)

#%%

exp_arr = cell_df_gene.to_numpy()

scaler = StandardScaler()
exp_arr = scaler.fit_transform(exp_arr)

pca = PCA(n_components=n_components)
exp_pcs = pca.fit_transform(exp_arr)

pc_mins = np.percentile(exp_pcs, min_pc_percentile, axis=0)
pc_maxs = np.percentile(exp_pcs, max_pc_percentile, axis=0)
exp_pcs = np.clip(exp_pcs, a_min=pc_mins, a_max=pc_maxs)

scaler = MinMaxScaler(feature_range=(0, 1))
exp_pcs = scaler.fit_transform(exp_pcs)

#%%

masks_im = tifffile.imread(masks_im_f)
masks_im = np.moveaxis(masks_im, 0, 2)

neur_im = tifffile.imread(neur_im_f)
neur_im = np.moveaxis(neur_im, 0, 2)

masks_im = filter_data.filt_masks_im_pp(masks_im, cell_df_orig)

#%%

pcs_im = np.ones(masks_im.shape + (n_components,), dtype=float)
pcs_im = pcs_im * 0.01

for i in range(n_components):
    pcs_im[:,:,:,i] = pcs_im[:,:,:,i] * neur_im

for i in tqdm(range(len(cell_df_orig))):
    
    z_plane = cell_df_orig['Z-plane'].iloc[i]
    pp = cell_df_orig['Principal Plane'].iloc[i]
    pp_python = z_plane + pp - 2
    
    mask = cell_df_orig['Mask ID'].iloc[i]
    mask_bool = (masks_im[:,:,pp_python] == mask)
    
    for j in range(n_components):
        pcs_im[mask_bool,pp_python,j] = neur_im[mask_bool,pp_python] * exp_pcs[i,j]
        
#%%

pcs_im = np.round(pcs_im).astype(np.uint16)
pcs_im = np.moveaxis(pcs_im, 2, 0)
pcs_im = np.moveaxis(pcs_im, 3, 1)
pcs_im = np.expand_dims(pcs_im, 0)

tifffile.imwrite(
    "pcs_rois.tif",
    pcs_im,
    imagej=True,         
    metadata={"axes": "TZCYX"},
)

#%%

anno_im_f = r'D:\2026-01-16_MC_SC_17\post\CCF_vols\s03L\s03Lannotation.tif'
anno_bound_im = img_utils.gen_anno_bound(anno_im_f, dilation_iter=3)
anno_bound_im = np.moveaxis(anno_bound_im, 2, 0)
anno_bound_im = np.expand_dims(anno_bound_im, 0)
anno_bound_im = anno_bound_im * 2000

# %%

tifffile.imwrite(
    "anno_bound.tif",
    anno_bound_im,
    imagej=True,         
    metadata={"axes": "TZYX"},
)



    
