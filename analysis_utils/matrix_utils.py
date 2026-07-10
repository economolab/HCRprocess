# -*- coding: utf-8 -*-
"""

"""

import bigfish.detection as detection
import numpy as np
import pandas as pd

import os

from analysis_utils import img_utils
from tqdm import tqdm

def derive_params(data_dir):
    
    (head,uniq_id) = os.path.split(data_dir)
    (post_dir,tail) = os.path.split(head)
    (head,tail) =  os.path.split(post_dir)
    (head,tail) = os.path.split(head)
    exp_name = tail[11:]
    
    return uniq_id, post_dir, exp_name

# Check if any tokens are in a string. Returns bool
def check_for_tokens(string, tokens):
    
    TF = any([token.casefold() in string.casefold() for token in tokens])
    
    return TF

def build_paths(data_dir, uniq_id, post_dir, exp_name, mode='curated'):
    
    reg_tokens = ['nt', 'neur', 'snap25']
    bin_tokens = ['bin']
    curated_mask_tokens = ['final_qc_masks']
    uncurated_mask_tokens = ['uncurated_masks']
    data_f = np.array(os.listdir(data_dir))
    
    reg_f_mask = [check_for_tokens(f,reg_tokens) for f in data_f]
    reg_f = data_f[reg_f_mask]
    data_f = data_f[np.invert(reg_f_mask)]
    
    bin_f_mask = [check_for_tokens(f,bin_tokens) for f in data_f]
    bin_f = data_f[bin_f_mask]
    if bin_f.size != 0:
        bin_f = str(bin_f[0])
    data_f = data_f[np.invert(bin_f_mask)]
    
    match mode:
        
        case 'curated':
            
            mask_f_mask = [check_for_tokens(f,curated_mask_tokens) for f in data_f]
            mask_f = data_f[mask_f_mask]
            mask_f = str(mask_f[0])
            data_f = data_f[np.invert(mask_f_mask)]
            
            mask_im_dir = os.path.join(post_dir,'masks',uniq_id)
            mask_im_dir_f = np.array(os.listdir(mask_im_dir))
            mask_im_f_mask = [check_for_tokens(f,['masks_qc_final']) for f in mask_im_dir_f]
            mask_im_f = mask_im_dir_f[mask_im_f_mask]
            mask_im_f = str(mask_im_f[0])
            mask_im_f = os.path.join(mask_im_dir, mask_im_f)
            
        case 'uncurated':
            
            mask_f_mask = [check_for_tokens(f,uncurated_mask_tokens) for f in data_f]
            mask_f = data_f[mask_f_mask]
            mask_f = str(mask_f[0])
            data_f = data_f[np.invert(mask_f_mask)]
            
            mask_im_dir = os.path.join(post_dir,'masks',uniq_id)
            mask_im_dir_f = np.array(os.listdir(mask_im_dir))
            mask_im_f_mask = [check_for_tokens(f,['cp_masks']) for f in mask_im_dir_f]
            mask_im_f = mask_im_dir_f[mask_im_f_mask]
            mask_im_f = str(mask_im_f[0])
            mask_im_f = os.path.join(mask_im_dir, mask_im_f)
            
    tifs_mask = [check_for_tokens(f,['.tif']) for f in data_f]
    data_f = data_f[tifs_mask]
    
    return data_f, reg_f, bin_f, mask_f, mask_im_f

def remove_label_channels(channels):
    
    channels = np.array(channels)
    labels = ['dtom', 'egfp', 'wga', 'ctb']
    
    TF = [check_for_tokens(channel, labels) for channel in channels]
    channels_filt = channels[~np.array(TF)]
    
    return channels_filt

def build_cell_df(mask_df, bin_df, remove_unassessed=True, remove_bad=True):
    
    cell_df = pd.DataFrame()
    cell_df.insert(0,'Cell ID',mask_df['Mask ID'])
    cell_df.insert(0,'Mask ID',mask_df['Mask ID'])
    cell_df.insert(cell_df.shape[1],'Z-plane',mask_df['Z-plane'])
    cell_df.insert(cell_df.shape[1],'Z-span',mask_df['Z-span'])
    cell_df.insert(cell_df.shape[1],'Principal Plane',mask_df['Principal Plane'])
    cell_df.insert(cell_df.shape[1],'Assessed',mask_df['Assessed'])
    cell_df.insert(cell_df.shape[1],'Result',mask_df['Result'])
    
    bin_df.drop(columns='Row',inplace=True)
    bin_genes = list(bin_df.columns)
    for bin_gene in bin_genes:
        cell_df.insert(cell_df.shape[1], bin_gene + '_b', bin_df[bin_gene])
        cell_df[bin_gene + '_b'] = cell_df[bin_gene + '_b'].astype('bool')
    
    if remove_unassessed:
        cell_df = cell_df[cell_df['Assessed'] == 1]
        
    if remove_bad:
        cell_df = cell_df[cell_df['Result'] == 'good']
    
    cell_df.drop(columns=['Assessed','Result'],inplace=True)
    cell_df.reset_index(drop=True,inplace=True)
    
    return cell_df

def count_spots(im, fluo, thresh=None):
    
    psf_sigma = img_utils.calc_psf_sigma(fluo)
    
    if thresh is not None:
        
        spots = detection.detect_spots(im, 
                                       threshold=thresh,
                                       voxel_size=325, 
                                       spot_radius=psf_sigma)
        
    else:
        
        spots = detection.detect_spots(im, 
                                       voxel_size=325, 
                                       spot_radius=psf_sigma)

    try:
        
        spots_post_decomposition, _, _ = detection.decompose_dense(
             im,
             spots=spots,
             voxel_size=325,
             spot_radius=psf_sigma)
        
    except RuntimeError:
        
        print('Unable to build reference spot: this plane is likely out of focus')
        spots_post_decomposition = spots
            
    return spots_post_decomposition


def find_zstack_thresh_bigfish(im, fluo):
    
    thresholds_bigfish = []
    psf_sigma = img_utils.calc_psf_sigma(fluo)
    
    for k in tqdm(range(im.shape[2])):
        _, thresh = detection.detect_spots(im[:,:,k], 
                                           return_threshold=True,
                                           voxel_size=325, 
                                           spot_radius=psf_sigma)
        thresholds_bigfish.append(thresh)
        
    thresholds_bigfish = np.array(thresholds_bigfish)
    
    return thresholds_bigfish
    
    
    
    


            
            