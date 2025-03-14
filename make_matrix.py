# -*- coding: utf-8 -*-
"""
Created on Fri Feb 28 14:02:44 2025

@author: jpv88
"""

import os

import numpy as np
import pandas as pd

from scipy.ndimage import distance_transform_edt
from skimage.filters import gaussian
from skimage.measure import regionprops
from skimage.restoration import rolling_ball
from PIL import Image, ImageSequence
from tqdm import tqdm

import bigfish.detection as detection
import bigfish.multistack as multistack
import bigfish.plot as plot

import matplotlib.pyplot as plt

from astropy.stats import SigmaClip
from photutils.background import Background2D, MedianBackground

# %% get all the folders and files and whatnot in order

label_image_path = r'D:\2024-11-19_MC-SC-4\post\masks\s01-left_cp_masks.tif'

head, tail = os.path.split(label_image_path)
uniq_id = tail.split('_')[0]
head, tail = os.path.split(head)

# intensity_image_folder = os.path.join(head, 'core_output', uniq_id)
intensity_image_folder = os.path.join(head, 'LipoGone', uniq_id)
intensity_image_files = os.listdir(intensity_image_folder)

reg_tokens = ["Snap25", "Snap", "NT"]

mask = np.ones(len(intensity_image_files), dtype=bool)
for i, intensity_image_file in enumerate(intensity_image_files):
    if any([el in intensity_image_file for el in reg_tokens]):
        mask[i] = False
        
intensity_image_files = np.array(intensity_image_files)[mask]
genes = [x.split('_')[1][:-7] for x in intensity_image_files]
genes_dict = {}
        
# %% load label image (Cellpose masks)

im = Image.open(label_image_path)

label_image = []
for page in ImageSequence.Iterator(im):
    label_image.append(np.array(page))
    
label_image = np.dstack(label_image)

cell_ids = list(np.unique(label_image))
cell_ids.remove(0)

cell_areas_dict = {}
for cell_id in cell_ids:
    cell_areas_dict[cell_id] = 0
    
unique, unique_counts = np.unique(label_image, return_counts=True)
unique = unique[1:]
unique_counts = unique_counts[1:]
for i, val in enumerate(unique):
    cell_areas_dict[val] = unique_counts[i]

# %% spot counting

# boolean mask that stores which cells will be kept (have no desaturated pixels)
preserve_label = np.ones(np.max(label_image), dtype=bool)

# iterate through genes
for i, intensity_image_file in tqdm(enumerate(intensity_image_files)):
    
    im = Image.open(os.path.join(intensity_image_folder, intensity_image_file))
    
    intensity_image = []
    for page in ImageSequence.Iterator(im):
        intensity_image.append(np.array(page))
        
    intensity_image = np.dstack(intensity_image)

    properties = regionprops(label_image, intensity_image=intensity_image)
    
    intensity_min = [x.intensity_min for x in properties]
    keep = (np.array(intensity_min) > 0)
    preserve_label = np.logical_and(preserve_label, keep)
    
    cell_spots_dict = {}
    for cell_id in cell_ids:
        cell_spots_dict[cell_id] = 0
        
    for j in tqdm(range(intensity_image.shape[2])):
        
        spots = detection.detect_spots(intensity_image[:,:,j], 
                                       voxel_size=(325, 325), 
                                       spot_radius=(325, 325))
        
        # spots_post_decomposition, dense_regions, reference_spot = detection.decompose_dense(
        #     intensity_image[:,:,j],
        #     spots=spots,
        #     voxel_size = (325, 325),
        #     spot_radius = (325, 325))
        
        spots_mask = np.zeros(np.shape(intensity_image[:,:,j]), dtype=bool)
        # for spot in spots_post_decomposition:
        for spot in spots:
            spots_mask[spot[0], spot[1]] = True
            
        label_spots = spots_mask * label_image[:,:,j]
        
        unique, unique_counts = np.unique(label_spots, return_counts=True)
        
        unique = unique[1:]
        unique_counts = unique_counts[1:]
        
        for k, val in enumerate(unique):
            cell_spots_dict[val] += unique_counts[k]
        
        
    genes_dict[genes[i]] = list(cell_spots_dict.values())
    
# %% mean intensity

# boolean mask that stores which cells will be kept (have no desaturated pixels)
preserve_label = np.ones(np.max(label_image), dtype=bool)

# iterate through genes
for i, intensity_image_file in tqdm(enumerate(intensity_image_files)):
    
    im = Image.open(os.path.join(intensity_image_folder, intensity_image_file))
    
    intensity_image = []
    for page in ImageSequence.Iterator(im):
        intensity_image.append(np.array(page))
        
    intensity_image = np.dstack(intensity_image)

    properties = regionprops(label_image, intensity_image=intensity_image)
    
    intensity_min = [x.intensity_min for x in properties]
    keep = (np.array(intensity_min) > 0)
    preserve_label = np.logical_and(preserve_label, keep)
    
    deltaFoverF = [x.intensity_mean for x in properties]

    background = np.zeros(intensity_image.shape, dtype=np.uint16)

    for j in tqdm(range(intensity_image.shape[2])):
        
        box_size = 50
        filter_size = 3

        sigma_clip = SigmaClip(sigma=3.0)
        bkg_estimator = MedianBackground()

        bkg = Background2D(intensity_image[:,:,j], box_size, filter_size=filter_size,
                            sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)
        
        background[:,:,j] = bkg.background

    bg_subtracted = np.subtract(intensity_image, background, dtype=np.int64)

    properties = regionprops(label_image, intensity_image=bg_subtracted)

    intensity_mean = [x.intensity_mean for x in properties]
    intensity_mean = np.array(intensity_mean)
    intensity_mean[intensity_mean < 0] = 0
    
    bg_properties = regionprops(label_image, intensity_image=background)
    
    bg_intensity_mean = [x.intensity_mean for x in bg_properties]
    bg_intensity_mean = np.array(bg_intensity_mean)
    
    deltaFoverF = intensity_mean / bg_intensity_mean
    
    genes_dict[genes[i]] = deltaFoverF

# %%

centroid = [x.centroid for x in properties]

x_cent = [x[0] for x in centroid]
y_cent = [x[1] for x in centroid]
z_cent = [x[2] for x in centroid]

x_cent = np.round(x_cent)
y_cent = np.round(y_cent)
z_cent = np.round(z_cent)

x_cent = x_cent.astype(np.uint16)
y_cent = y_cent.astype(np.uint16)
z_cent = z_cent.astype(np.uint16)

genes_dict["x_cent"] = x_cent
genes_dict["y_cent"] = y_cent
genes_dict["z_cent"] = z_cent

genes_dict["area"] = list(cell_areas_dict.values())
genes_dict["cell_id"] = cell_ids
    
# %%

gene_df = pd.DataFrame(genes_dict)
gene_df = gene_df[preserve_label]

gene_df.to_csv(uniq_id + "_df.csv", index=False)

# %%

# coords = [x.coords for x in properties]

# bg = np.zeros(len(coords))
# for j, coords_label in tqdm(enumerate(coords)):
    
#     bg_coords = []
#     for coord in coords_label:
#         bg_coords.append(indices[:,coord[0],coord[1],coord[2]])
#     bg_coords = set(map(tuple,bg_coords))
    
#     bg_label = []
#     for bg_coord in bg_coords:
#         bg_label.append(intensity_image[bg_coord[0],bg_coord[1],bg_coord[2]])
        
#     bg[j] = np.median(bg_label)

# filtered_intensity_image = np.zeros(intensity_image.shape, dtype=float)

# # iterate through z-planes
# for j in tqdm(range(intensity_image.shape[2])):
#     filtered_intensity_image[:,:,j] = gaussian(intensity_image[:,:,j], 
#                                                sigma=3,
#                                                preserve_range=True)
    
# filtered_intensity_image = np.round(filtered_intensity_image)
# filtered_intensity_image = filtered_intensity_image.astype(np.uint16)

# background = np.zeros(intensity_image.shape, dtype=np.uint16)

# for j in tqdm(range(filtered_intensity_image.shape[2])):
    
#     box_size = 50
#     filter_size = 3

#     sigma_clip = SigmaClip(sigma=3.0)
#     bkg_estimator = MedianBackground()

#     bkg = Background2D(filtered_intensity_image[:,:,j], box_size, filter_size=filter_size,
#                        sigma_clip=sigma_clip, bkg_estimator=bkg_estimator)
    
#     background[:,:,j] = bkg.background

# bg_subtracted = np.subtract(intensity_image, background, dtype=np.int64)

# properties = regionprops(label_image, intensity_image=bg_subtracted)

# intensity_mean = [x.intensity_mean for x in properties]
# intensity_mean = np.array(intensity_mean)
# intensity_mean[intensity_mean < 0] = 0


# for i in range(len(spots_z)):
#     z_idx = np.repeat(i, len(spots_z[i]))
#     spots_test[i] = np.hstack((np.expand_dims(z_idx, 1), spots_z[i]))
    
# for j in range(len(spots_z)):
#     z_idx = np.repeat(j, len(spots_z[j]))
#     spots_z[j] = np.hstack((np.expand_dims(z_idx, 1), spots_z[j]))

# filtered_intensity_image = np.zeros(intensity_image.shape, dtype=float)

# iterate through z-planes
# for j in tqdm(range(intensity_image.shape[2])):
#     filtered_intensity_image[:,:,j] = gaussian(intensity_image[:,:,j], 
#                                                 sigma=3,
#                                                 preserve_range=True)
    
# filtered_intensity_image = np.round(filtered_intensity_image)
# filtered_intensity_image = filtered_intensity_image.astype(np.uint16)

# %%



test_im = intensity_image[:,:,20]
test_im_label = label_image[:,:,20]

spots = detection.detect_spots(test_im, 
                               voxel_size=(325, 325), 
                               spot_radius=(325, 325))

spots_post_decomposition, dense_regions, reference_spot = detection.decompose_dense(
    test_im,
    spots=spots,
    voxel_size = (325, 325),
    spot_radius = (325, 325))

fov_results = multistack.extract_cell(test_im_label, 2, rna_coord=spots_post_decomposition)
df = multistack.summarize_extraction_results(fov_results, 2)

plot.plot_detection(test_im, spots_post_decomposition, contrast=True)
plot.plot_reference_spot(reference_spot, rescale=True)

# %%

norm_rna = df["nb_rna"] / df["cell_area"]

# %%

plt.scatter(gene_df["area"], gene_df["Ebf3"])
# plt.scatter(gene_df["z_cent"], gene_df["Ebf3"])



