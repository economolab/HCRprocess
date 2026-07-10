# -*- coding: utf-8 -*-
"""

"""

import tifffile

import matplotlib.pyplot as plt
import numpy as np

from bigfish.stack import log_filter
from scipy import ndimage
from sklearn.linear_model import RANSACRegressor
from skimage.morphology import disk
from skimage.segmentation import find_boundaries
from tqdm import tqdm

def gen_anno_bound(anno_im_f, dilation_iter=1):
    
    anno_im = tifffile.imread(anno_im_f)
    anno_im = np.moveaxis(anno_im, 0, 2)
    
    anno_bound_im = np.zeros(np.shape(anno_im), dtype=np.uint16)
    for i in tqdm(range(np.shape(anno_im)[2])):
        anno_bound_im[:,:,i] = find_boundaries(anno_im[:,:,i])
        
    for i in tqdm(range(np.shape(anno_bound_im)[2])):
        anno_bound_im[:,:,i] = ndimage.binary_dilation(anno_bound_im[:,:,i], 
                                                       structure=disk(1),
                                                       iterations=dilation_iter)
    
    return anno_bound_im

# source: Gaussian approximations of fluorescence microscope point-spread function models
# Bo Zhang, Josiane Zerubia, and Jean-Christophe Olivo-Marin (2007)
# PSF sigma when modeling as a Gaussian (in same units as lambda_em, should be nm)
def calc_psf_sigma(fluo, NA=1.15):
    
    fluo_to_lambda_em = {
        "Atto 425": 484,
        "Alexa Fluor 488": 525,
        "Alexa Fluor 546": 572,
        "Alexa Fluor 594": 618,
        "Alexa Fluor 647": 671
    }
    
    lambda_em = fluo_to_lambda_em[fluo]
    sigma = 0.21 * (lambda_em / NA)
    
    return sigma

def correct_depth_attenuation(im, plot_fit=False, scale_to_log_filter=False, fluo=None, log_filter_zstack=False):
    
    if scale_to_log_filter and fluo is None:
        raise ValueError("Must pass fluorophore argument when scaling correction to log filter")
        
    im = np.float64(im)
    z_indices = np.arange(im.shape[2])
    z_medians = np.median(im, axis=(0, 1))
    
    log_medians = np.zeros_like(z_medians)
    nonzero_mask = (z_medians != 0)
    log_medians[nonzero_mask] = np.log(z_medians[nonzero_mask])
    
    reg = RANSACRegressor(min_samples=0.5, max_trials=int(1e4))
    reg.fit(z_indices.reshape(-1, 1), log_medians)
    alpha = -reg.estimator_.coef_
    
    correction = np.exp(alpha * z_indices)
    corrected_im = im * correction[np.newaxis, np.newaxis, :]
    
    if plot_fit:
        
        plt.figure()
        plt.xlabel("Z-index")
        plt.ylabel("log of median intensity")
        plt.plot(z_indices, log_medians, label='data')
        plt.scatter(z_indices[~reg.inlier_mask_], log_medians[~reg.inlier_mask_], c='r', marker='x')
        plt.scatter(z_indices[reg.inlier_mask_], log_medians[reg.inlier_mask_], c='g', marker='o')
        plt.plot(z_indices, reg.predict(z_indices.reshape(-1, 1)), label='fitted')
        plt.legend()
        
        corrected_z_means = np.mean(corrected_im, axis=(0, 1))
        
        plt.figure()
        plt.xlabel("Z-index")
        plt.ylabel("Median intensity")
        plt.plot(z_indices, z_medians, label='data')
        plt.plot(z_indices, corrected_z_means, label='corrected')
        plt.legend()
        
    if scale_to_log_filter:
        
        psf_sigma = calc_psf_sigma(fluo)
        voxel_size = 325
        sigma = psf_sigma / voxel_size
        
        n_planes = corrected_im.shape[2]
        im_filtered = np.zeros_like(corrected_im)
        
        uint16_max = np.iinfo(np.uint16).max
        
        restart = True
        while restart:
            
            restart = False
            
            match log_filter_zstack:
                
                case False:
                    
                    for k in tqdm(range(n_planes)):
                        
                        im_filtered[:,:,k] = log_filter(corrected_im[:,:,k], sigma)
                        plane_max = np.max(im_filtered[:,:,k])
                        ceiling_broken = (plane_max > uint16_max)
                        
                        if ceiling_broken:
                            print(f"Uint16 ceiling broken by value: {plane_max}, decrementing 1%")
                            corrected_im *= 0.99
                            restart = True
                            break
                        
                case True:
                    
                    im_filtered = log_filter(corrected_im, sigma)
                    stack_max = np.max(im_filtered)
                    ceiling_broken = (stack_max > uint16_max)
                    
                    if ceiling_broken:
                        print(f"Uint16 ceiling broken by value: {stack_max}, decrementing 1%")
                        corrected_im *= 0.99
                        restart = True
                        break
                    
    corrected_im = np.round(corrected_im).astype(np.uint16)
        
    return corrected_im

def make_crop_slices(pp_masks_im, margin=100):
    
    pp_masks_im_b = (pp_masks_im > 0)
    pp_masks_im_b = pp_masks_im_b.astype(np.uint8)
    pp_masks_im_b_max = np.max(pp_masks_im_b, axis=2)
    im = pp_masks_im_b_max
    
    im_n_rows = im.shape[0]
    im_n_cols = im.shape[1]
    
    slices = ndimage.find_objects(im)[0]
    
    row_slice = slice(max(slices[0].start - margin, 0), min(slices[0].stop + margin, im_n_rows), None)
    col_slice = slice(max(slices[1].start - margin, 0), min(slices[1].stop + margin, im_n_cols), None)
    margin_slices = (row_slice, col_slice)
    
    return margin_slices

