# -*- coding: utf-8 -*-
"""

"""

#%%

import tifffile

import numpy as np

from sklearn.preprocessing import StandardScaler, MinMaxScaler
from sklearn.decomposition import PCA
from tqdm import tqdm

from analysis_utils import load_data, filter_data, img_utils, matrix_utils

#%%

exp_dirs = [r'D:\2026-01-16_MC_SC_17']
# masks_im_f = r'D:\2026-01-16_MC_SC_17\post\masks\s03L\s03L_masks_qc_final.tif'
masks_im_f = r'D:\2026-01-16_MC_SC_17\post\masks\s03L\s03L_cp_masks.tif'
neur_im_f = r'D:\2026-01-16_MC_SC_17\post\core_output\s03L\r1_neur445.tif'

n_components = 3
min_pc_percentile = 10
max_pc_percentile = 90

from scipy import ndimage
    
def make_projection_image(masks_im, neur_im, cell_df, latents, 
                          background_intensity=0.1):
    
    n_components = latents.shape[1]
    
    print('Allocating memory for projection image...')
    proj_im = neur_im[..., np.newaxis].astype(float)
    proj_im = proj_im * background_intensity
    proj_im = np.repeat(proj_im, n_components, axis=-1)
    
    z_planes = cell_df['Z-plane'].to_numpy()
    pps = cell_df['Principal Plane'].to_numpy()
    pp_python_arr = (z_planes + pps - 2).astype(int)
    mask_ids = cell_df['Mask ID'].to_numpy()
    
    print('Writing planes in projection image...')
    for z in tqdm(np.unique(pp_python_arr)):
        mask_plane = masks_im[:, :, z]
        neur_plane = neur_im[:, :, z]

        objects = ndimage.find_objects(mask_plane)
        bbox_by_label = {label: bbox for label, bbox in enumerate(objects, start=1) if bbox is not None}

        rows = np.where(pp_python_arr == z)[0]
        for i in rows:
            
            mask_id = int(mask_ids[i])
            
            if mask_id not in bbox_by_label:
                raise ValueError(
                    f"mask_id={mask_id} not found in z-plane {z}'s mask image "
                    f"(cell row {i})"
                )
                
            bbox = bbox_by_label[mask_id]

            sub_mask = (mask_plane[bbox] == mask_id)
            sub_neur = neur_plane[bbox]

            y_idx, x_idx = np.nonzero(sub_mask)         
            y_idx = y_idx + bbox[0].start                  
            x_idx = x_idx + bbox[1].start
        
            proj_im[y_idx, x_idx, z, :] = sub_neur[sub_mask][:, None] * latents[i, :]
        
    return proj_im

def write_projection_image(proj_im_f, proj_im):
    
    print('Saving tiff...')
    proj_im = np.round(proj_im).astype(np.uint16)
    proj_im = np.moveaxis(proj_im, 2, 0) # move z to the front
    proj_im = np.moveaxis(proj_im, 3, 1) # move c to the second dimension
    proj_im = np.expand_dims(proj_im, 0) # add an empty t dimension at the front

    tifffile.imwrite(
        proj_im_f,
        proj_im,
        bigtiff=True,    
        ome=True,
        metadata={"axes": "TZCYX"},
    )

    
#%%

# cell_df_orig = load_data.load_experiments(exp_dirs)
# cell_df_orig = filter_data.filt_cells(cell_df_orig, 's03L')

# cell_df_gene = load_data.load_experiments(exp_dirs)
# cell_df_gene = filter_data.filt_cells(cell_df_gene, 's03L')
# cell_df_gene = filter_data.filt_features(cell_df_gene, remove_features=['dtom', 'egfp'], genes=True)

# cell_df_orig = load_data.load_cell_df('MC_SC_17_s03L_uncurated.csv')
# cell_df_orig = filter_data.filt_cells(cell_df_orig, 's03L')

# cell_df_gene = load_data.load_cell_df('MC_SC_17_s03L_uncurated.csv')
# cell_df_gene = filter_data.filt_cells(cell_df_gene, 's03L')
# cell_df_gene = filter_data.filt_features(cell_df_gene, genes=True)

cell_df_orig = load_data.load_cell_df('MC_SC_17_s03L_new_spot_all_cells.csv')
cell_df_orig = filter_data.filt_cells(cell_df_orig, 's03L')

cell_df_gene = load_data.load_cell_df('MC_SC_17_s03L_new_spot_all_cells.csv')
cell_df_gene = filter_data.filt_cells(cell_df_gene, 's03L')
cell_df_gene = filter_data.filt_features(cell_df_gene, genes=True)
cell_df_gene = filter_data.filt_features(cell_df_gene, remove_features=['dtom_b', 'egfp_b'], genes=True)

#%%

masks_im = tifffile.imread(masks_im_f)
masks_im = np.moveaxis(masks_im, 0, 2)

neur_im = tifffile.imread(neur_im_f)
neur_im = np.moveaxis(neur_im, 0, 2)

masks_im = filter_data.filt_masks_im_pp(masks_im, cell_df_orig)

#%%

latents = matrix_utils.embed_cell_df(cell_df_gene, embed='umap', n_components=3, ref_mask=None)
proj_im = make_projection_image(masks_im, neur_im, cell_df_orig, latents)
write_projection_image('pca_projection.tif', proj_im)


#%%

nonzero_cells = (np.sum(cell_df_gene.to_numpy(), axis=1) != 0)
good_masks = np.uint16(np.loadtxt('IRNPARNmasks.csv'))
good_masks_bool = np.isin(cell_df_orig['Mask ID'].values, good_masks)

cell_df_gene = cell_df_gene[nonzero_cells]
cell_df_orig = cell_df_orig[nonzero_cells]
good_masks_bool = good_masks_bool[nonzero_cells]

#%%

latents = matrix_utils.embed_cell_df(cell_df_gene, embed='umap', n_components=3, ref_mask=good_masks_bool)
proj_im = make_projection_image(masks_im, neur_im, cell_df_orig, latents)
write_projection_image('umap_projection.tif', proj_im)

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

proj_im = make_projection_image(masks_im, neur_im, cell_df_orig, exp_pcs)
write_projection_image('pca_projection.tif', proj_im)


#%%

pcs_im = np.ones(masks_im.shape + (n_components,), dtype=float)
pcs_im = pcs_im * 0.1

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
pcs_im = np.moveaxis(pcs_im, 2, 0) # move z to the front
pcs_im = np.moveaxis(pcs_im, 3, 1)
pcs_im = np.expand_dims(pcs_im, 0)

tifffile.imwrite(
    "pcs_rois.tif",
    pcs_im,
    imagej=True,         
    metadata={"axes": "TZCYX"},
)

#%%

nonzero_cells = (np.sum(cell_df_gene.to_numpy(), axis=1) != 0)
good_masks = np.uint16(np.loadtxt('IRNPARNmasks.csv'))
good_masks_bool = np.isin(cell_df_orig['Mask ID'].values, good_masks)

cell_df_gene = cell_df_gene[nonzero_cells]
cell_df_orig = cell_df_orig[nonzero_cells]
good_masks_bool = good_masks_bool[nonzero_cells]

#%%

latents = matrix_utils.embed_cell_df(cell_df_gene, embed='pca', n_components=3, ref_mask=good_masks_bool)
proj_im = make_projection_image(masks_im, neur_im, cell_df_orig, latents)
write_projection_image('pca_projection.tif', proj_im)

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
#%%
def make_projection_image(masks_im, neur_im, cell_df, latents, 
                          background_intensity=0.1):
    
    n_components = latents.shape[1]
    proj_im = np.ones(masks_im.shape + (n_components,), dtype=float)
    proj_im = proj_im * background_intensity
    
    for i in range(n_components):
        proj_im[:,:,:,i] = proj_im[:,:,:,i] * neur_im
    
    for i in tqdm(range(len(cell_df))):
        
        z_plane = cell_df['Z-plane'].iloc[i]
        pp = cell_df['Principal Plane'].iloc[i]
        pp_python = z_plane + pp - 2
        
        mask = cell_df['Mask ID'].iloc[i]
        mask_bool = (masks_im[:,:,pp_python] == mask)
        
        for j in range(n_components):
            proj_im[mask_bool,pp_python,j] = neur_im[mask_bool,pp_python] * latents[i,j]
        
    return proj_im



    
