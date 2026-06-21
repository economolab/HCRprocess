# -*- coding: utf-8 -*-
"""
Take a fully completed HCRprocess image directory and turn it into a data
matrix.
"""

# %% imports

import os
import tifffile

import numpy as np
import pandas as pd

from tqdm import tqdm

from analysis_utils import matrix_utils, filter_data, img_utils

#%% explicit and derived params

data_dir = 'D:\\2026-01-16_MC_SC_17\\post\\core_output\\s03L'

uniq_id, post_dir, exp_name = matrix_utils.derive_params(data_dir)
data_f, reg_f, bin_f, mask_f, mask_im_f = matrix_utils.build_paths(data_dir, 
                                                                   uniq_id, 
                                                                   post_dir, 
                                                                   exp_name,
                                                                   mode='curated')

#%% load everything

mask_df = pd.read_csv(os.path.join(data_dir,mask_f))
bin_df = pd.read_csv(os.path.join(data_dir,bin_f))
cell_df = matrix_utils.build_cell_df(mask_df, bin_df)

data_f = matrix_utils.remove_label_channels(data_f)

quant_genes = [f[3:-7] for f in data_f]
quant_genes.remove('slc17a6')
quant_genes.remove('slc32a1')
quant_genes.remove('slc5a7')
data_f = data_f[data_f != 'r1_slc17a6488.tif']
data_f = data_f[data_f != 'r2_slc32a1488.tif']
data_f = data_f[data_f != 'r4_slc5a7488.tif']

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
    
masks_im = tifffile.imread(mask_im_f)
masks_im = np.transpose(masks_im,[1,2,0])

pp_masks_im = filter_data.filt_masks_im_pp(masks_im, cell_df)
row_slice, col_slice = img_utils.make_crop_slices(pp_masks_im)
pp_masks_im_crop = pp_masks_im [row_slice, col_slice, :]

gene_to_fluo = {
    "phox2b": "Alexa Fluor 594",
    "ralyl": "Alexa Fluor 647",
    "slc17a6": "Alexa Fluor 488",
    "tenm2": "Alexa Fluor 546",
    "ebf3": "Alexa Fluor 647",
    "pcp4": "Alexa Fluor 594",
    "slc32a1": "Alexa Fluor 488",
    "tshz2": "Alexa Fluor 546",
    "alcam": "Alexa Fluor 546",
    "celf2": "Alexa Fluor 647",
    "meis2": "Alexa Fluor 594",
    "rph3a": "Alexa Fluor 488",
    "robo1": "Alexa Fluor 647",
    "slc5a7": "Alexa Fluor 488",
    "syt1": "Alexa Fluor 594",
    "zfhx3": "Alexa Fluor 546"
    }

alcam_thresh = np.load("alcam_thresh_s03L.npy")
alcam_thresh = 6660
celf2_thresh = np.load("celf2_thresh_s03L.npy")
celf2_thresh = 4400
ebf3_thresh = np.load("ebf3_thresh_s03L.npy")
ebf3_thresh = 4840
meis2_thresh = np.load("meis2_thresh_s03L.npy")
meis2_thresh = 5540
pcp4_thresh = np.load("pcp4_thresh_s03L.npy")
pcp4_thresh = 6400
phox2b_thresh = np.load("phox2b_thresh_s03L.npy")
phox2b_thresh = 6330
ralyl_thresh = np.load("ralyl_thresh_s03L.npy")
ralyl_thresh = 4420
robo1_thresh = np.load("robo1_thresh_s03L.npy")
robo1_thresh = 3915
rph3a_thresh = np.load("rph3a_thresh_s03L.npy")
rph3a_thresh = 7160
syt1_thresh = np.load("syt1_thresh_s03L.npy")
syt1_thresh = 5870
tenm2_thresh = np.load("tenm2_thresh_s03L.npy")
tenm2_thresh = 4064
tshz2_thresh = np.load("tshz2_thresh_s03L.npy")
tshz2_thresh = 5290
zfhx3_thresh = np.load("zfhx3_thresh_s03L.npy")
zfhx3_thresh = 5640

thresholds = [phox2b_thresh, ralyl_thresh, tenm2_thresh, ebf3_thresh, pcp4_thresh,
              tshz2_thresh, alcam_thresh, celf2_thresh, meis2_thresh, rph3a_thresh,
              robo1_thresh, syt1_thresh, zfhx3_thresh]

def plot_thresholds(thresholds, gene):
    
    plt.plot(thresholds)
    plt.xlabel("Z-index")
    plt.ylabel("Bigfish threshold")
    plt.title(gene)

#%% threshold finding

i = 6
quant_gene = 'alcam'

fluo = gene_to_fluo[quant_gene]
intensity_image = np.copy(quant_gene_ims[i])
intensity_image_crop = intensity_image[row_slice, col_slice, :]
intensity_image_crop = img_utils.correct_depth_attenuation(intensity_image_crop, 
                                                      scale_to_log_filter=True,
                                                      plot_fit=True,
                                                      fluo=fluo)

# thresholds_bigfish = matrix_utils.find_zstack_thresh_bigfish(intensity_image_crop, fluo)


# %% spot counting

import bigfish.detection as detection

# iterate through genes
for i, quant_gene in enumerate(quant_genes):
    
    print('Quantifying ' + quant_gene + '...')
    
    fluo = gene_to_fluo[quant_gene]
    intensity_image = np.copy(quant_gene_ims[i])
    intensity_image_crop = intensity_image[row_slice, col_slice, :]
    intensity_image_crop = img_utils.correct_depth_attenuation(intensity_image_crop, 
                                                          scale_to_log_filter=True, 
                                                          fluo=fluo)
    
    cell_spots_dict = {}
    for cell_id in cell_df['Mask ID'].values:
        cell_spots_dict[cell_id] = 0
        
    for k in tqdm(range(intensity_image_crop.shape[2])):
        
        spots = matrix_utils.count_spots(intensity_image_crop[:,:,k], fluo, thresh=thresholds[i])
        
        spots_mask = np.zeros(np.shape(intensity_image_crop[:,:,k]), dtype=bool)
        
        if spots is not None:
            for spot in spots:
                spots_mask[spot[0], spot[1]] = True
            
        label_spots = spots_mask * pp_masks_im_crop[:,:,k]
        
        unique, unique_counts = np.unique(label_spots, return_counts=True)
        unique = unique[1:]
        unique_counts = unique_counts[1:]
        
        if len(unique) != 0:
            for idx, val in enumerate(unique):
                cell_spots_dict[val] += unique_counts[idx]
        
    cell_df[quant_gene] = list(cell_spots_dict.values())
    
#%%

cell_df.to_csv(exp_name + '_' + uniq_id + '_new_spot' + '.csv',
               index=False)
    
#%%

import matplotlib.pyplot as plt
import bigfish.detection as detection

from bigfish.plot import plot_detection
from bigfish.stack import log_filter
from scipy import ndimage

import pwlf

from scipy.ndimage import uniform_filter1d
from scipy.stats import median_abs_deviation

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

def set_init_guess(im, fluo):
    
    spot_radius = calc_psf_sigma(fluo)
    
    spots, threshold = detection.detect_spots(im, 
                                   return_threshold=True,
                                   voxel_size=(325, 325), 
                                   spot_radius=spot_radius)
    
    return threshold

def find_min_thresh(im, x0, fluo, size_init_range, max_spots=1e5, step=100):
    
    spot_radius = calc_psf_sigma(fluo)
    current_thresh = x0
    num_spots = 0
    
    with tqdm(desc='Finding minimum threshold') as pbar:
        
        while num_spots < max_spots:
            
            current_thresh -= step
            
            if current_thresh < 1:
                min_thresh = 1
                return min_thresh
            
            if current_thresh < (x0 - size_init_range/2):
                min_thresh = (x0 - size_init_range / 2)
                return min_thresh
            
            spots = detection.detect_spots(im, 
                                           threshold=current_thresh,
                                           voxel_size=(325, 325), 
                                           spot_radius=spot_radius)
            
            num_spots = spots.shape[0]
        
            pbar.set_postfix(current_thresh=f"{current_thresh:.1f}", num_spots=num_spots)
            pbar.update(1)
        
    min_thresh = current_thresh
    
    return min_thresh

def set_init_range(im, fluo, size_init_range=2e3):
    
    x0 = set_init_guess(im, fluo)
    min_thresh = find_min_thresh(im, x0, fluo, size_init_range)
    init_range = np.linspace(min_thresh, min_thresh + size_init_range, 100, dtype=float)
    
    return init_range

def set_range(x0, min_thresh, size_range):
    
    bounds = [x0 - size_range/2, x0 + size_range/2]
    
    if bounds[0] < min_thresh:
        bounds[0] = min_thresh
        bounds[1] = min_thresh + size_range
    
    search_range = np.linspace(bounds[0], bounds[1], 100)
    
    return search_range 

def build_counts_curve(im, thresh_range, fluo):
    
    spot_radius = calc_psf_sigma(fluo)
    counts = []

    for thresh in tqdm(thresh_range):
        
        spots = detection.detect_spots(im, 
                                       threshold=thresh,
                                       voxel_size=(325, 325), 
                                       spot_radius=spot_radius)
    
        count = spots.shape[0]
        counts.append(count)
        
    return counts

def preprocess_counts(counts):
    
    counts = [np.float64(i) for i in counts]
    deriv1 = np.gradient(counts)
    deriv1 = uniform_filter1d(deriv1, 5)
    deriv1 = np.abs(deriv1)
    
    return deriv1

def calc_fano(windows):
    
    num_windows = windows.shape[0]
    fano = np.zeros(num_windows)
    
    zero_windows = np.all((windows == 0), axis=1)
    fano[~zero_windows] = np.var(windows[~zero_windows], axis=1) / np.mean(windows[~zero_windows], axis=1)
    
    return fano

def calc_med_sugg(fano, thresholds_fano):
    
    median = np.median(fano)
    MAD = median_abs_deviation(fano)
    
    min_target = median - MAD
    max_target = median + MAD
    
    min_sugg_idx = np.argmax(fano <= min_target)
    max_sugg_idx = np.argmax(fano <= max_target)
    
    med_suggs = [thresholds_fano[min_sugg_idx], thresholds_fano[max_sugg_idx]]
    med_sugg = np.mean(med_suggs)
    
    return med_sugg

def calc_pwlf_sugg(fano, thresholds_fano):
    
    my_pwlf = pwlf.PiecewiseLinFit(thresholds_fano, fano)
    res = my_pwlf.fit(2)
    
    intercept_sugg = res[1]
    
    x1 = res[1]
    x2 = res[2]
    y1 = my_pwlf.predict(x1)
    y2 = my_pwlf.predict(x2)
    slope = (y2 - y1) / (x2 - x1)
    intercept = y1 - (slope * x1)
    fit_fano = (slope * thresholds_fano) + intercept
    
    diff = fano - fit_fano
    intersections = np.where(np.diff(np.sign(diff)))[0]
    first_intersect = intersections[0] if len(intersections) > 0 else None
    
    if first_intersect is not None:
        intersect_sugg = thresholds_fano[first_intersect]
    else:
        intersect_sugg = intercept_sugg
    
    if intersect_sugg < intercept_sugg:
        intersect_sugg = intercept_sugg
        
    return intercept_sugg, intersect_sugg

def suggest_thresh_window(window_shape, deriv1, thresholds):
    
    windows = np.lib.stride_tricks.sliding_window_view(deriv1, window_shape)
    idx = int((window_shape - 1) / 2)
    thresholds_fano = thresholds[idx:-idx]
        
    fano = calc_fano(windows)
    med_sugg = calc_med_sugg(fano, thresholds_fano)
    intercept_sugg, intersect_sugg = calc_pwlf_sugg(fano, thresholds_fano)
    threshold_sugg = np.mean((med_sugg, intercept_sugg, intersect_sugg))
    
    return threshold_sugg

def suggest_thresh(counts, thresholds):
    
    deriv1 = preprocess_counts(counts)
    window_shapes = [3, 5, 7, 9, 11, 13, 15, 17, 19, 21]

    threshold_suggs = []
    for window_shape in window_shapes:
        threshold_suggs.append(suggest_thresh_window(window_shape, deriv1, thresholds))

    threshold_sugg = np.mean(threshold_suggs)
    
    return threshold_sugg

# source: Gaussian approximations of fluorescence microscope point-spread function models
# Bo Zhang, Josiane Zerubia, and Jean-Christophe Olivo-Marin (2007)
# PSF sigma when modeling as a Gaussian (in same units as lambda_em, should be nm)
def calc_psf_sigma(fluo, NA=1.15):
    
    fluo_to_lambda_em = {
        "Alexa Fluor 546": 572,
        "Alexa Fluor 594": 618,
        "Alexa Fluor 647": 671
    }
    
    lambda_em = fluo_to_lambda_em[fluo]
    sigma = 0.21 * (lambda_em / NA)
    
    return sigma

def find_thresh(im, fluo):
    
    init_range = set_init_range(im, fluo)
    min_thresh = init_range[0]
    
    counts = build_counts_curve(im, init_range, fluo)
    threshold_sugg = suggest_thresh(counts, init_range)
    print(threshold_sugg)
    
    # range_sizes = [500, 200, 100, 50, 20]
    
    # for range_size in range_sizes:
        
    #     search_range = set_range(threshold_sugg, min_thresh, range_size)
    #     counts = build_counts_curve(im, search_range)
    #     threshold_sugg = suggest_thresh(counts, search_range)
    #     print(threshold_sugg)
    
    thresh = round(threshold_sugg)
    
    return thresh

def count_spots(im, thresh, fluo):
    
    spot_radius = calc_psf_sigma(fluo)
    
    spots = detection.detect_spots(im, 
                                   threshold=thresh,
                                   voxel_size=(325, 325), 
                                   spot_radius=spot_radius)
        
    spots_post_decomposition, dense_regions, reference_spot = detection.decompose_dense(
         im,
         spots=spots,
         voxel_size = (325, 325),
         spot_radius = spot_radius)
    
    return spots_post_decomposition

def count_spots_zstack(im, slices, threshes, fluo):
    
    spot_radius = calc_psf_sigma(fluo)
    spots_z = []
    
    for k in tqdm(range(im.shape[2])):
        
        im_k = im[:,:,k]
        im_k_crop = im_k[slices]
    
        spots = detection.detect_spots(im_k_crop, 
                                       threshold=threshes[k],
                                       voxel_size=(325, 325), 
                                       spot_radius=spot_radius)
        
        try:
            spots_post_decomposition, dense_regions, reference_spot = detection.decompose_dense(
                 im_k_crop,
                 spots=spots,
                 voxel_size = (325, 325),
                 spot_radius = spot_radius)
        except RuntimeError:
            print(f'Unable to build reference spot: plane {k} is likely out of focus')
            spots_post_decomposition = None
                
        spots_z.append(spots_post_decomposition)
    
    return spots_z

def filter_spots(spots, masks_im):
    
    masks_im_b = (masks_im > 0)
    spots_filt = spots[masks_im_b[spots[:,0],spots[:,1]],:]
    
    return spots_filt

def find_thresh_gene(im, slices=None):
    
    num_planes = im.shape[2]
    
    thresholds = []
    for k in range(num_planes):
        
        im_k = im[:,:,k]
        if slices is not None:
            im_k = im_k[slices]
        
        thresh = find_thresh(im_k)
        thresholds.append(thresh)
        
    return thresholds



#%%

slices = make_crop_slices(pp_masks_im)

meis2_thresh = find_thresh_gene(quant_gene_ims[10], slices=slices)
ralyl_thresh = find_thresh_gene(quant_gene_ims[1], slices=slices)
zfhx3_thresh = find_thresh_gene(quant_gene_ims[15], slices=slices)

# %%

thresholds = []
for j in tqdm(range(10, 21)):
    intensity_image_j = quant_gene_ims[10][:,:,j]
    intensity_image_j_crop = intensity_image_j[slices]
    thresholds.append(find_thresh(intensity_image_j_crop, "Alexa Fluor 594"))
        

# %%

slices = make_crop_slices(pp_masks_im)

j = 20
intensity_image_j = quant_gene_ims[10][:,:,j]
intensity_image_j_crop = intensity_image_j[slices]

min_thresh = find_min_thresh(intensity_image_j_crop, 5126, "Alexa Fluor 594")
search_range = set_range(5126, 4126, 2000)
counts = build_counts_curve(intensity_image_j_crop, search_range, "Alexa Fluor 594")

# %%

thresholds_bigfish = []
for j in tqdm(range(quant_gene_ims[10].shape[2])):
    intensity_image_j = quant_gene_ims[10][:,:,j]
    intensity_image_j_crop = intensity_image_j[slices]
    thresholds_bigfish.append(set_init_guess(intensity_image_j_crop, "Alexa Fluor 594"))
    
thresholds_bigfish = np.array(thresholds_bigfish)
    
#%%

def oned_filter(data):
    
    filt_data = uniform_filter1d(data, 3)
    
    return filt_data

# def find_outliers(residuals, thresh=2):
    
#     std = np.std(residuals)
#     outliers = residuals > thresh*std
    
#     return outliers

def purge_outliers(x, y, y_pred, thresh=0.2):
    
    residuals = (np.abs(y_pred - y) / y_pred)
    outliers = residuals > thresh
    
    new_x = x[~outliers]
    new_y = y[~outliers]
    
    return new_x, new_y, outliers
    
def fit_poly(x, y, degree=4):
    
    p = np.polyfit(x, y, degree)
    fit_val = np.polyval(p, x)
    
    return fit_val

def find_zstack_thresh_bigfish(im, slices, fluo):
    
    thresholds_bigfish = []
    for j in tqdm(range(im.shape[2])):
        intensity_image_j = im[:,:,j]
        intensity_image_j_crop = intensity_image_j[slices]
        thresholds_bigfish.append(set_init_guess(intensity_image_j_crop, fluo))
        
    thresholds_bigfish = np.array(thresholds_bigfish)
    
    old_z = np.array(list(range(im.shape[2])))
    old_thresh = thresholds_bigfish
    
    while True:
        
        fit_val = fit_poly(old_z, old_thresh)
        
        plt.figure()
        plt.plot(old_z, old_thresh)
        plt.plot(old_z, fit_val)
        
        new_z, new_thresh, outliers = purge_outliers(old_z, old_thresh, fit_val)

        print(np.sum(outliers))        
        if np.sum(outliers) == 0:
            break
        
        old_z = new_z
        old_thresh = new_thresh

    p = np.polyfit(new_z, new_thresh, 4)
    z = np.array(list(range(im.shape[2])))
    final_thresh = np.polyval(p, z)
    
    return thresholds_bigfish, final_thresh

thresholds_bigfish, final_thresh = find_zstack_thresh_bigfish(quant_gene_ims[1], slices, "Alexa Fluor 647")

#%%

from scipy.optimize import curve_fit

def find_zstack_thresh_bigfish(im, slices, fluo):
    
    thresholds_bigfish = []
    for j in tqdm(range(im.shape[2])):
        intensity_image_j = im[:,:,j]
        intensity_image_j_crop = intensity_image_j[slices]
        thresholds_bigfish.append(set_init_guess(intensity_image_j_crop, fluo))
        
    thresholds_bigfish = np.array(thresholds_bigfish)
    
    return thresholds_bigfish

def exp_decay(x, a, b, c):
    return a * np.exp(-b * x) + c

def fix_bigfish_thresh_simple(thresholds_bigfish):
    
    n_planes = len(thresholds_bigfish)
    z_indices = np.arange(n_planes)
    
    mean_thresh = np.mean(thresholds_bigfish)
    # mask = (thresholds_bigfish > mean_thresh)
    mask = [4, 6, 8, 13, 16, 18, 21, 28, 32]
    
    thresholds_bigfish_filt = thresholds_bigfish[mask]
    z_indices_filt = z_indices[mask]
    
    p0 = (3000, 0.2, 4000)
    popt, pcov = curve_fit(exp_decay, z_indices_filt, thresholds_bigfish_filt, p0=p0)
    
    final_thresh = exp_decay(z_indices, popt[0], popt[1], popt[2])
    
    return final_thresh
    
    
    

# %%

spots_z = count_spots_zstack(quant_gene_ims[1], slices, final_thresh, "Alexa Fluor 647")

#%%

from skimage.exposure import match_histograms

def histogram_match_stack(stack, reference_slice_idx=0):
    
    reference = stack[:,:,reference_slice_idx]
    corrected = np.zeros_like(stack)
    
    for z in tqdm(range(stack.shape[2])):
        corrected[:,:,z] = match_histograms(stack[:,:,z], reference)
        
    return corrected

test = histogram_match_stack(quant_gene_ims[1], reference_slice_idx=7)

# %%
thresholds_bigfish, final_thresh = find_zstack_thresh_bigfish(test, slices, "Alexa Fluor 647")
spots_z = count_spots_zstack(test, slices, final_thresh, "Alexa Fluor 647")

#%%

from sklearn.base import BaseEstimator, RegressorMixin
from sklearn.linear_model import RANSACRegressor

# class ExponentialDecayRegressor(BaseEstimator, RegressorMixin):
        
#     def fit(self, X, y):
        
#         X = np.squeeze(X)
#         self.coeffs_ = np.polyfit(X, y, 1)
        
#         return self
        
#     def predict(self, X):
        
#         X = np.squeeze(X)
#         y = np.polyval(self.coeffs_, X)
        
#         return y
    
def correct_depth_attenuation(im, plot_fit=False):
    
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
    corrected_im = np.uint16(corrected_im)
    
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
        plt.ylabel("Mean intensity")
        plt.plot(z_indices, z_medians, label='data')
        plt.plot(z_indices, corrected_z_means, label='corrected')
        plt.legend()
        
    return corrected_im

from bigfish.stack import log_filter
from bigfish.detection import local_maximum_detection, spots_thresholding

def log_filter_zstack(im, fluo):
    
    psf_sigma = calc_psf_sigma(fluo)
    voxel_size = 325
    sigma = psf_sigma / voxel_size
    
    n_planes = im.shape[2]
    im_filtered = np.zeros_like(im)
    
    for k in tqdm(range(n_planes)):
        im_filtered[:,:,k] = log_filter(im[:,:,k], sigma)
        
    return im_filtered

def local_maximum_detection_zstack(im, fluo):
    
    psf_sigma = calc_psf_sigma(fluo)
    voxel_size = 325
    min_distance = psf_sigma / voxel_size
    
    n_planes = im.shape[2]
    mask = np.zeros_like(im, dtype=bool)
    
    for k in tqdm(range(n_planes)):
        mask[:,:,k] = local_maximum_detection(im[:,:,k], min_distance)
        
    return mask

def get_thresh_range_zstack(im, num_thresh=100):
    
    start_ranges = np.percentile(im_filtered, 95, axis=(0, 1))
    start_range = np.median(start_ranges)
    end_ranges = np.percentile(im_filtered, 99.9999, axis=(0, 1))
    end_range = np.median(end_ranges)
    thresh_range = np.linspace(start_range, end_range, num=num_thresh)
    
    return thresh_range

def spots_thresholding_zstack(im, mask, thresh_range):
    
    n_planes = im.shape[2]
    num_thresh = len(thresh_range)
    spots_per_thresh = np.zeros((n_planes, num_thresh))
    
    for k in tqdm(range(n_planes)):
        
        j = 0
        
        for thresh in tqdm(thresh_range):
            
            spots, _ = spots_thresholding(im[:,:,k], mask[:,:,k], thresh)
            spots_per_thresh[k,j] = spots.shape[0]
            j += 1
            
    return spots_per_thresh

def detect_spots_zstack(im, mask, final_thresh, fluo):
    
    n_planes = im.shape[2]
    spots_per_plane = []
    
    spot_radius = calc_psf_sigma(fluo)
    
    for k in tqdm(range(n_planes)):
        
        spots, _ = spots_thresholding(im[:,:,k], mask[:,:,k], final_thresh[k])
        
        try:
            spots_post_decomposition, dense_regions, reference_spot = detection.decompose_dense(
                 im[:,:,k],
                 spots=spots,
                 voxel_size = (325, 325),
                 spot_radius = spot_radius)
        except RuntimeError:
            print(f'Unable to build reference spot: plane {k} is likely out of focus')
            spots_post_decomposition = None
        
        spots_per_plane.append(spots_post_decomposition)

    return spots_per_plane

def get_thresh(thresholds, spot_counts):
    
    my_pwlf = pwlf.PiecewiseLinFit(thresholds, spot_counts)
    res = my_pwlf.fit(2)
    
    thresh = res[1]
    
    return thresh

def count_spots_zstack(im, threshes, fluo):
    
    spot_radius = calc_psf_sigma(fluo)
    spots_z = []
    
    for k in tqdm(range(im.shape[2])):
        
        im_k = im[:,:,k]
    
        spots = detection.detect_spots(im_k, 
                                       threshold=threshes[k],
                                       voxel_size=(325, 325), 
                                       spot_radius=spot_radius)
        
        try:
            spots_post_decomposition, dense_regions, reference_spot = detection.decompose_dense(
                 im_k,
                 spots=spots,
                 voxel_size = (325, 325),
                 spot_radius = spot_radius)
        except RuntimeError:
            print(f'Unable to build reference spot: plane {k} is likely out of focus')
            spots_post_decomposition = None
                
        spots_z.append(spots_post_decomposition)
    
    return spots_z


        

    
    
        
# %%
slices = make_crop_slices(pp_masks_im)
im = im[slices[0], slices[1], :]
corrected_im = correct_depth_attenuation(im, plot_fit=True)

#%%

import pandas as pd

fluo = "Alexa Fluor 594"

im_filtered = log_filter_zstack(corrected_im, fluo)
mask = local_maximum_detection_zstack(im_filtered, fluo)
thresh_range = get_thresh_range_zstack(im_filtered)
spots_per_thresh = spots_thresholding_zstack(im_filtered, mask, thresh_range)

final_thresh = []
for i in range(spots_per_thresh.shape[0]):
    final_thresh.append(get_thresh(thresh_range, spots_per_thresh[i,:]))
    
plt.plot(final_thresh)
final_thresh = pd.Series(final_thresh)
final_thresh = final_thresh.rolling(3, center=True, min_periods=1)
final_thresh = final_thresh.mean().to_numpy()
plt.plot(final_thresh)

spots_z = detect_spots_zstack(im_filtered, mask, final_thresh, fluo)
plt.plot(spots_z)

#%%
plt.figure()
plt.imshow(im_filtered[:,:,13], vmax=20000)
plt.figure()
plt.imshow(im_filtered[:,:,14], vmax=20000)
#%%

spots_per_thresh = []
for thresh in tqdm(thresh_range):
    spots, _ = spots_thresholding(im_filtered[:,:,20], mask[:,:,20], thresh)
    spots_per_thresh.append(spots.shape[0])
    


#%%

data = spots_per_thresh
fig, ax = plt.subplots(figsize=(10, 5))

colors = plt.cm.plasma(np.linspace(0, 1, data.shape[0]))
for row, color in zip(data, colors):
    ax.plot(row, color=color)

sm = plt.cm.ScalarMappable(cmap="plasma", norm=plt.Normalize(0, data.shape[0] - 1))
fig.colorbar(sm, ax=ax, label="Z plane")

ax.set_xlabel("Threshold")
ax.set_ylabel("Value")
plt.tight_layout()
plt.show()

#%%

from scipy.interpolate import UnivariateSpline
from sklearn.isotonic import IsotonicRegression
from scipy.signal import savgol_filter


def curvature(x, y):
    
    dy = np.gradient(y)
    
    iso_reg = IsotonicRegression(increasing=True)
    iso_reg.fit(x, dy)
    dy_mono = iso_reg.predict(x)
    dy_savgol = savgol_filter(dy_mono, 35, 3)
    spl = UnivariateSpline(thresh_range, dy_savgol, s=0.01)
    
    d2y = spl.derivative(1)(thresh_range)
    
    kappa = np.abs(d2y) / (1 + spl(thresh_range**2)**1.5)
    
    return kappa

dy = np.gradient(spots_per_thresh[20,:])

iso_reg = IsotonicRegression(increasing=True)
iso_reg.fit(thresh_range, dy)
dy_mono = iso_reg.predict(thresh_range)
dy_savgol = savgol_filter(dy_mono, 35, 3)
spl = UnivariateSpline(thresh_range, dy_savgol, s=0.01)


plt.plot(dy)
plt.plot(dy_mono)
plt.plot(thresh_range, dy_savgol)
plt.plot(thresh_range, np.arctan(dy_savgol))
plt.plot(spl.derivative(1)(thresh_range))

d2y = np.gradient(dy_savgol)
plt.plot(d2y)


#%%

from kneed import KneeLocator

kl = KneeLocator(thresh_range, 
                 spots_per_thresh[20,:], 
                 curve='concave',
                 direction='decreasing')


#%%
        
        
    

# im = quant_gene_ims[1]

# z_means = np.mean(im, axis=(0, 1))
# z_indices = np.arange(im.shape[2])
# z_indices = z_indices[7:]

# log_means = np.log(z_means)
# log_means = log_means[7:]

# coeffs = np.polyfit(z_indices, log_means, 1)
# alpha = coeffs[0]

# z_indices = np.arange(im.shape[2])
# correction = np.exp(alpha * z_indices)

# corrected_im = im / correction[np.newaxis, np.newaxis, :]
# corrected_im = corrected_im.astype(np.uint16)

# %%
thresholds_bigfish = find_zstack_thresh_bigfish(corrected_im, slices, "Alexa Fluor 647")
spots_z = count_spots_zstack(corrected_im, slices, final_thresh, "Alexa Fluor 647")

# %%

def count_spots_per_plane(spots_z):
    
    num_spots = []
    for spots in spots_z:
        if spots is None:
            num_spots.append(0)
        else:
            num_spots.append(spots.shape[0])
    
    return num_spots

num_spots = count_spots_per_plane(spots_z)

#%%

ninetynineth_percentiles = []
for k in tqdm(range(corrected_im.shape[2])):
    sigma = calc_psf_sigma("Alexa Fluor 647")
    filt_plane = log_filter(corrected_im[:,:,k],sigma/325)
    ninetynineth_percentiles.append(np.percentile(filt_plane, 99.9999))
#%%
    
thresh = set_init_guess(intensity_image_j_crop, "Alexa Fluor 594")
test = find_min_thresh(intensity_image_j_crop, thresh)

# xHat = np.linspace(min(thresholds_fano), max(thresholds_fano), num=10000)
# yHat = my_pwlf.predict(xHat)
# plt.figure()
# plt.plot(thresholds_fano, fano, '-')
# plt.plot(xHat, yHat, '-')
# plt.show()


#%%

spots = count_spots(intensity_image_j_crop, 4671, "Alexa Fluor 594")

pp_masks_im_j = pp_masks_im[:,:,j]
pp_masks_im_j_crop = pp_masks_im_j[slices]

spots = filter_spots(spots, pp_masks_im_j_crop)

plot_detection(intensity_image_j_crop, spots, 
               rescale=True, contrast=False, radius=1)


#%%


    
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
            
        label_spots = spots_mask * pp_masks_im[:,:,j]
        
        unique, unique_counts = np.unique(label_spots, return_counts=True)
        
        unique = unique[1:]
        unique_counts = unique_counts[1:]
        
        for k, val in enumerate(unique):
            cell_spots_dict[val] += unique_counts[k]
        
        
    cell_df[quant_gene] = list(cell_spots_dict.values())
    
#%%

tifffile.imwrite(exp_name + '_' + uniq_id + '_pp_masks.tif', 
                 np.moveaxis(pp_masks_im,2,0), 
                 photometric='minisblack')

#%%

# cell_df.to_csv(exp_name + '_' + uniq_id + '.csv',
#                index=False)

cell_df.to_csv(exp_name + '_' + uniq_id + '_uncurated' + '.csv',
               index=False)

# #%% quantify genes 

# for i, quant_gene in enumerate(quant_genes):
                          
#     print('Quantifying ' + quant_gene + '...')
#     intensity_image = np.copy(quant_gene_ims[i])
#     background = np.zeros(intensity_image.shape, dtype=np.uint16)
    
#     for j in tqdm(range(intensity_image.shape[2]),
#                   desc='Calculating background for ' + quant_gene):
        
#         box_size = 50
#         filter_size = 3
    
#         sigma_clip = SigmaClip(sigma=3.0)
#         bkg_estimator = MedianBackground()
    
#         bkg = Background2D(intensity_image[:,:,j], box_size, filter_size=filter_size,
#                             sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)
        
#         background[:,:,j] = bkg.background
    
    
#     intensity_image = (intensity_image - background) / background
#     props = regionprops(pp_mask_im,intensity_image=intensity_image)
    
#     intensity_mean = [prop.intensity_mean for prop in props]
#     cell_df[quant_gene] = intensity_mean
    
# #%% quantify genes, no bg

# for i, quant_gene in enumerate(quant_genes):
                          
#     print('Quantifying ' + quant_gene + '...')
#     intensity_image = np.copy(quant_gene_ims[i])
    
#     props = regionprops(pp_mask_im,intensity_image=intensity_image)
    
#     intensity_mean = [prop.intensity_mean for prop in props]
#     cell_df[quant_gene] = intensity_mean

# %%

# slices = make_crop_slices(pp_masks_im)

# from skimage import filters

# def tenengrad_xy(im):
    
#     filt_im = filters.sobel(im)
#     tenengrad = np.mean(np.square(filt_im)) / np.mean(im)
    
#     return tenengrad

# def tenengrad_z(im, slices=None):
    
#     n_planes = im.shape[2]
#     tenengrad = []
    
#     for k in tqdm(range(n_planes)):
        
#         if slices is not None:
#             input_im = im[:,:,k][slices]
#         else:
#             input_im = im[:,:,k]
            
#         tenengrad.append(tenengrad_xy(input_im))
    
#     return tenengrad

# tenengrad = tenengrad_z(quant_gene_ims[10], slices=slices)

# # %%


# def laplacian_var_xy(im):
    
#     filt_im = ndimage.laplace(im)
#     laplace_var = np.var(filt_im)
    
#     return laplace_var

# def laplacian_var_z(im, slices=None):
    
#     n_planes = im.shape[2]
#     laplace_var = []
    
#     for k in tqdm(range(n_planes)):
        
#         if slices is not None:
#             input_im = im[:,:,k][slices]
#         else:
#             input_im = im[:,:,k]
            
#         laplace_var.append(laplacian_var_xy(input_im))
    
#     return laplace_var

# laplace_var = laplacian_var_z(quant_gene_ims[1], slices=slices)

# def find_good_planes(im, fluo, slices=None):
    
#     spot_radius = calc_psf_sigma(fluo)
    
#     for k in tqdm(range(im.shape[2])):
        
#         im_k = im[:,:,k]
        
#         if slices is not None:
#             im_k = im_k[slices]
    
#         spots = detection.detect_spots(im_k,
#                                        threshold=4000,
#                                        voxel_size=(325, 325), 
#                                        spot_radius=spot_radius)
            
#         spots_post_decomposition, dense_regions, reference_spot = detection.decompose_dense(
#              im_k,
#              spots=spots,
#              voxel_size = (325, 325),
#              spot_radius = spot_radius)
