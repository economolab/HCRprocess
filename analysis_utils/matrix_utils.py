# -*- coding: utf-8 -*-
"""

"""

import numpy as np

import os

def derive_params(data_dir):
    
    (head,uniq_id) = os.path.split(data_dir)
    (post_dir,tail) = os.path.split(head)
    (head,tail) =  os.path.split(post_dir)
    (head,tail) = os.path.split(head)
    exp_name = tail[11:]
    
    return uniq_id, post_dir, exp_name

# Check if any tokens are in a string. Returns bool
def check_for_tokens(string,tokens):
    
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


            
            