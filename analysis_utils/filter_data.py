# -*- coding: utf-8 -*-
"""

"""

import numpy as np
import pandas as pd

from tqdm import tqdm

# filter a cell df to only cells whose cell ids contain a target keyword
def filt_cells(cell_df, kw):
    
    contains_kw = [kw in cell_id for cell_id in cell_df.index]
    cell_df = cell_df[contains_kw]
    
    return cell_df

# filter features from a cell df
def filt_features(cell_df, remove_features=None, keep_features=None, genes=False):
    
    if genes == True:
        
        non_gene_cols = ['Mask ID', 'Z-plane', 'Z-span', 'Principal Plane', 'Area']
        remove_cols = [col in non_gene_cols for col in cell_df.columns]
        keep_cols = [not x for x in remove_cols]
        cell_df = cell_df.loc[:,keep_cols]
    
    if remove_features != None:
        
        remove_cols = [col in remove_features for col in cell_df.columns]
        keep_cols = [not x for x in remove_cols]
        cell_df = cell_df.loc[:,keep_cols]
        
    if keep_features != None:
        
        keep_cols = [col in keep_features for col in cell_df.columns]
        cell_df = cell_df.loc[:,keep_cols]
    
    return cell_df

# Filter a masks image to only contain principal planes of masks in a 
# corresponding cell dataframe
def filt_masks_im_pp(masks_im, cell_df):
    
    cell_df['Mask ID'] = cell_df['Mask ID'].astype('uint16')
    cell_df['Z-plane'] = cell_df['Z-plane'].astype('uint16')
    cell_df['Principal Plane'] = cell_df['Principal Plane'].astype('uint16')
    
    # all the masks in the masks image
    all_masks = np.unique(masks_im)
    all_masks = np.delete(all_masks,0)
    
    # which masks in the masks image are also in the cell dataframe, and 
    # should therefore be kept
    good_masks = cell_df['Mask ID'].values
    good_masks_bool = [mask in good_masks for mask in all_masks]
    
    # which masks in the masks image are not in the cell dataframe, and
    # should therefore be cut
    masks_to_delete_bool = np.invert(good_masks_bool)
    masks_to_delete = all_masks[masks_to_delete_bool]
    
    # deleting the bad masks
    isin_bool = np.isin(masks_im, masks_to_delete)
    masks_im[isin_bool] = 0

    pp_masks_im = np.zeros(masks_im.shape, dtype=np.uint16)
    
    for mask in tqdm(good_masks):
        
        z_plane = cell_df['Z-plane'][cell_df['Mask ID'] == mask].values[0]
        pp = cell_df['Principal Plane'][cell_df['Mask ID'] == mask].values[0]
        pp_python = z_plane + pp - 2
        pp_bool = (masks_im[:,:,pp_python] == mask)
        
        pp_masks_im[pp_bool,pp_python] = mask
        
    return pp_masks_im