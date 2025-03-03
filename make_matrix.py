# -*- coding: utf-8 -*-
"""
Created on Fri Feb 28 14:02:44 2025

@author: jpv88
"""

import os

import numpy as np
import pandas as pd

from scipy.ndimage import distance_transform_edt
from skimage.filters import threshold_otsu
from skimage.measure import regionprops
from PIL import Image, ImageSequence
from tqdm import tqdm

# %% get all the folders and files and whatnot in order

label_image_path = r'D:\2024-11-19_MC-SC-4\post\masks\s01-left_cp_masks.tif'

head, tail = os.path.split(label_image_path)
uniq_id = tail.split('_')[0]
head, tail = os.path.split(head)

intensity_image_folder = os.path.join(head, 'core_output', uniq_id)
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

bin_label_image = (label_image > 0)
indices = distance_transform_edt(bin_label_image,
                                 return_distances=False,
                                 return_indices=True)

# %% 

preserve_label = np.ones(np.max(label_image), dtype=bool)

for i, intensity_image_file in tqdm(enumerate(intensity_image_files)):
    
    im = Image.open(os.path.join(intensity_image_folder, intensity_image_file))
    
    intensity_image = []
    for page in ImageSequence.Iterator(im):
        intensity_image.append(np.array(page))
        
    intensity_image = np.dstack(intensity_image)

    properties = regionprops(label_image, intensity_image=intensity_image)
    
    intensity_mean = [x.intensity_mean for x in properties]
    
    intensity_min = [x.intensity_min for x in properties]
    keep = (np.array(intensity_min) > 0)
    preserve_label = np.logical_and(preserve_label, keep)
    
    coords = [x.coords for x in properties]
    
    bg = np.zeros(len(coords))
    for j, coords_label in tqdm(enumerate(coords)):
        
        bg_coords = []
        for coord in coords_label:
            bg_coords.append(indices[:,coord[0],coord[1],coord[2]])
        bg_coords = set(map(tuple,bg_coords))
        
        bg_label = []
        for bg_coord in bg_coords:
            bg_label.append(intensity_image[bg_coord[0],bg_coord[1],bg_coord[2]])
            
        bg[j] = np.median(bg_label)
    
    genes_dict[genes[i]] = (intensity_mean - bg)
    
    # genes_dict[genes[i]] = intensity_mean
    
# %%

gene_df = pd.DataFrame(genes_dict)
gene_df = gene_df[preserve_label]

gene_df.to_csv(uniq_id + "_df.csv")

# %%


