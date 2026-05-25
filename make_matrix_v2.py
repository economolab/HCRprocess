# -*- coding: utf-8 -*-
"""
Take a fully completed HCRprocess image directory and turn it into a data
matrix.
"""

# %% imports

import os

import numpy as np
import pandas as pd

import tifffile

from astropy.stats import SigmaClip
from photutils.background import Background2D, MedianBackground
from skimage.measure import regionprops
from tqdm import tqdm

from analysis_utils import matrix_utils

#%% explicit and derived params

data_dir = 'D:\\2026-01-16_MC_SC_17\\post\\core_output\\s03L'

uniq_id, post_dir, exp_name = matrix_utils.derive_params(data_dir)
data_f, reg_f, bin_f, mask_f, mask_im_f = matrix_utils.build_paths(data_dir, 
                                                                   uniq_id, 
                                                                   post_dir, 
                                                                   exp_name,
                                                                   mode='uncurated')

#%% load everything

# bin_df = pd.read_csv(os.path.join(data_dir,bin_f))
# bin_df.drop(columns='Row',inplace=True)
# bin_genes = list(bin_df.columns)
# bin_genes.remove('rph3a')

# bin_genes_f_mask = [check_for_tokens(f,bin_genes) for f in data_f]
# bin_genes_f = data_f[bin_genes_f_mask]
# data_f = data_f[np.invert(bin_genes_f_mask)]

data_f = data_f[2:]

mask_df = pd.read_csv(os.path.join(data_dir,mask_f))

cell_df = pd.DataFrame()
cell_df.insert(0,'Cell ID',mask_df['Mask ID'])
cell_df.insert(0,'Mask ID',mask_df['Mask ID'])
cell_df.insert(cell_df.shape[1],'Z-plane',mask_df['Z-plane'])
cell_df.insert(cell_df.shape[1],'Z-span',mask_df['Z-span'])
cell_df.insert(cell_df.shape[1],'Principal Plane',mask_df['Principal Plane'])
cell_df.insert(cell_df.shape[1],'Assessed',mask_df['Assessed'])
cell_df.insert(cell_df.shape[1],'Result',mask_df['Result'])

# for bin_gene in bin_genes:
#     cell_df.insert(cell_df.shape[1],bin_gene,bin_df[bin_gene])
#     cell_df[bin_gene] = cell_df[bin_gene].astype('bool')
    
cell_df = cell_df[cell_df['Assessed'] == 1]
cell_df = cell_df[cell_df['Result'] == 'good']
cell_df.drop(columns=['Assessed','Result'],inplace=True)
cell_df.reset_index(drop=True,inplace=True)

quant_genes = [f[3:-7] for f in data_f]

for quant_gene in quant_genes:
    cell_df.insert(cell_df.shape[1],quant_gene,np.zeros(cell_df.shape[0]))
    
cell_df['Cell ID'] = cell_df['Cell ID'].astype(str)
for i in range(cell_df.shape[0]):
    cell_df.loc[i,"Cell ID"] = exp_name + '_' + uniq_id + '_' + str(i+1)

    
#%% load images

quant_gene_ims = []

for i in tqdm(range(len(data_f)),desc='Loading images...'):
    im = tifffile.imread(os.path.join(data_dir,data_f[i]))
    im = np.transpose(im,[1,2,0])
    im[np.isnan(im)] = 0
    quant_gene_ims.append(im)
    
mask_im = tifffile.imread(mask_im_f)
mask_im = np.transpose(mask_im,[1,2,0])
    
#%%

# delete masks from masks image that didn't make the cut
all_masks = np.unique(mask_im)
all_masks = np.delete(all_masks,0)
good_masks = cell_df['Mask ID'].values
good_masks_bool = [mask in good_masks for mask in all_masks]
masks_to_delete_bool = np.invert(good_masks_bool)
masks_to_delete = all_masks[masks_to_delete_bool]
isin_bool = np.isin(mask_im,masks_to_delete)
mask_im[isin_bool] = 0

pp_mask_im = np.zeros(mask_im.shape, dtype=np.uint16)
for mask in tqdm(good_masks):
    z_plane = cell_df['Z-plane'][cell_df['Mask ID'] == mask].values[0]
    pp = cell_df['Principal Plane'][cell_df['Mask ID'] == mask].values[0]
    pp_python = z_plane + pp - 2
    pp_bool = (mask_im[:,:,pp_python] == mask)
    pp_mask_im[pp_bool,pp_python] = mask
    
props = regionprops(pp_mask_im)
area = [prop.area for prop in props]
cell_df.insert(cell_df.shape[1],'Area',area)

#%% quantify genes 

for i, quant_gene in enumerate(quant_genes):
                          
    print('Quantifying ' + quant_gene + '...')
    intensity_image = np.copy(quant_gene_ims[i])
    background = np.zeros(intensity_image.shape, dtype=np.uint16)
    
    for j in tqdm(range(intensity_image.shape[2]),
                  desc='Calculating background for ' + quant_gene):
        
        box_size = 50
        filter_size = 3
    
        sigma_clip = SigmaClip(sigma=3.0)
        bkg_estimator = MedianBackground()
    
        bkg = Background2D(intensity_image[:,:,j], box_size, filter_size=filter_size,
                            sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)
        
        background[:,:,j] = bkg.background
    
    
    intensity_image = (intensity_image - background) / background
    props = regionprops(pp_mask_im,intensity_image=intensity_image)
    
    intensity_mean = [prop.intensity_mean for prop in props]
    cell_df[quant_gene] = intensity_mean
    
#%% quantify genes, no bg

for i, quant_gene in enumerate(quant_genes):
                          
    print('Quantifying ' + quant_gene + '...')
    intensity_image = np.copy(quant_gene_ims[i])
    
    props = regionprops(pp_mask_im,intensity_image=intensity_image)
    
    intensity_mean = [prop.intensity_mean for prop in props]
    cell_df[quant_gene] = intensity_mean
    
# %% spot counting

import bigfish.detection as detection

# iterate through genes
for i, quant_gene in enumerate(quant_genes):
    
    print('Quantifying ' + quant_gene + '...')
    intensity_image = np.copy(quant_gene_ims[i])
    
    cell_spots_dict = {}
    for cell_id in cell_df['Mask ID'].values:
        cell_spots_dict[cell_id] = 0
        
    for j in tqdm(range(intensity_image.shape[2])):
        
        spots = detection.detect_spots(intensity_image[:,:,j], 
                                       voxel_size=(325, 325), 
                                       spot_radius=(325, 325))
        
        spots_post_decomposition, dense_regions, reference_spot = detection.decompose_dense(
            intensity_image[:,:,j],
            spots=spots,
            voxel_size = (325, 325),
            spot_radius = (325, 325))
        
        spots_mask = np.zeros(np.shape(intensity_image[:,:,j]), dtype=bool)
        for spot in spots_post_decomposition:
        # for spot in spots:
            spots_mask[spot[0], spot[1]] = True
            
        label_spots = spots_mask * pp_mask_im[:,:,j]
        
        unique, unique_counts = np.unique(label_spots, return_counts=True)
        
        unique = unique[1:]
        unique_counts = unique_counts[1:]
        
        for k, val in enumerate(unique):
            cell_spots_dict[val] += unique_counts[k]
        
        
    cell_df[quant_gene] = list(cell_spots_dict.values())
    
#%%

tifffile.imwrite(exp_name + '_' + uniq_id + '_pp_masks.tif', 
                 np.moveaxis(pp_mask_im,2,0), 
                 photometric='minisblack')

#%%

# cell_df.to_csv(exp_name + '_' + uniq_id + '.csv',
#                index=False)

cell_df.to_csv(exp_name + '_' + uniq_id + '_uncurated' + '.csv',
               index=False)
