# -*- coding: utf-8 -*-
"""
Shading correction for images collected on spinning disk confocal
"""

import nd2
import os
import tifffile

import matplotlib.pyplot as plt
import numpy as np

from basicpy import BaSiC
from pathlib import Path
from tqdm import tqdm

import sys

from scipy.optimize import curve_fit
from scipy.ndimage import uniform_filter1d

#%%

# path = sys.argv[1]
path = r'D:\2026-07-13_HCR_WGAc_YH95\2026-07-16_r1_HCR\raw\s01L__slc17a6488_ralyl647_phox2b594_tenm2546_neur445__HCR.nd2'

print('Loading image file...')
myfile = nd2.ND2File(path)
my_array = nd2.imread(path)

p = Path(path)
parent = str(p.parent)
stem = str(p.stem)

#%% collect metadata

experiment = myfile.experiment

for loop in experiment:
    if loop.type == 'XYPosLoop':
        XYPosLoop = loop
        
XYPosLoopParams = XYPosLoop.parameters
points = XYPosLoopParams.points

names = [point.name for point in points]
stagePositionsUm = [point.stagePositionUm for point in points]
x = [stagePositionUm.x for stagePositionUm in stagePositionsUm]
y = [stagePositionUm.y for stagePositionUm in stagePositionsUm]

# naming convention is (letter, number) / (row, index)
rows = [name[0] for name in names]
cols = [name[1] for name in names]

# 1-indexed
rows = [ord(row) - ord('A') + 1 for row in rows]
cols = [int(col) for col in cols]

n_rows = max(rows)
n_cols = max(cols)
n_pos = myfile.sizes['P']
n_planes = myfile.sizes['Z']
n_chans = myfile.sizes['C']

#%% 

# photobleaching types
# type 1: no bleach
# type 2: left bleach only
# type 3: top bleach only
# type 4: top bleach and right bleach
# type 5: top bleach and left bleach

bleach_types = []

for i in range(n_pos):
    
    # very first tile
    if rows[i] == 1 and cols[i] == 1:
        bleach_types.append(1)
    
    # first row but not the first tile
    elif rows[i] == 1:
        bleach_types.append(2)
        
    # even row, last column
    elif rows[i] % 2 == 0 and cols[i] == n_cols:
        bleach_types.append(3)
    
    # odd row, first column
    elif rows[i] % 2 != 0 and cols[i] == 1:
        bleach_types.append(3)
        
    # even row
    elif rows[i] % 2 == 0:
        bleach_types.append(4)
    
    # odd row
    elif rows[i] % 2 != 0:
        bleach_types.append(5)
        
mask_1 = (np.array(bleach_types) == 1)
mask_2 = (np.array(bleach_types) == 2)
mask_3 = (np.array(bleach_types) == 3)
mask_4 = (np.array(bleach_types) == 4)
mask_5 = (np.array(bleach_types) == 5)

indices_1 = np.arange(n_pos)[mask_1]
indices_2 = np.arange(n_pos)[mask_2]
indices_3 = np.arange(n_pos)[mask_3]
indices_4 = np.arange(n_pos)[mask_4]
indices_5 = np.arange(n_pos)[mask_5]

left_bleach_mask = mask_2 | mask_5
indices_left_bleach = np.arange(n_pos)[left_bleach_mask]

right_bleach_mask = mask_4
indices_right_bleach = np.arange(n_pos)[right_bleach_mask]

top_bleach_mask = mask_3 | mask_4 | mask_5
indices_top_bleach = np.arange(n_pos)[top_bleach_mask]

#%%

# fit basic object for a channel
def fit_channel_basic(image, channel_idx, smoothness_flatfield=5):
    
    basic_obj = BaSiC(get_darkfield=False, 
                      fitting_mode='ladmap', 
                      sort_intensity=True,
                      smoothness_flatfield=smoothness_flatfield)
    
    image_chan = image[:,:,channel_idx,:,:]
    t, z, y, x = image_chan.shape
    image_chan = image_chan.reshape(t*z, y, x)
    
    basic_obj.fit(image_chan)
        
    return basic_obj

def correct_channel_basic(image, basic_obj, channel_idx):
  
    for i in tqdm(range(n_pos)):
        for j in range(n_planes):
            input_image = image[i,j,channel_idx,:,:]
            input_image = np.expand_dims(input_image, axis=0)
            image[i,j,channel_idx,:,:] = basic_obj.transform(input_image, use_tqdm=False)
        
    return image

def normalize_by_median(img):
    
    if isinstance(img, list):
        img_norm = [im / np.median(im) for im in img]
    else:
        img_norm = img / np.median(img)
        
    return img_norm

def compute_median_flatfield(images):
    
    images = normalize_by_median(images)
    stacked_images = np.stack(images, axis=0)
    flatfield = np.median(stacked_images, axis=0)
    
    return flatfield

def compute_intensity_profiles(img, debug_plot=False):
    
    n_rows, n_cols = img.shape
    
    x_profile = np.median(img, axis=0)
    y_profile = np.median(img, axis=1)
    
    x_profile = uniform_filter1d(x_profile, round(0.01*n_cols), mode='nearest')
    y_profile = uniform_filter1d(y_profile, round(0.01*n_rows), mode='nearest')
    
    if debug_plot:
        
        fig, axes = plt.subplots(1, 3)
        
        axes[0].imshow(img)
        axes[1].plot(x_profile)
        axes[2].plot(y_profile)
        
        axes[0].set_title("Flat-field image")
        axes[1].set_title("X profile (left to right)")
        axes[2].set_title("Y profile (top to bottom)")
        
        axes[0].set_box_aspect(1)
        axes[1].set_box_aspect(1)
        axes[2].set_box_aspect(1)
    
    return x_profile, y_profile

def gaussian(x, a, b, c):
    
    y = a * np.exp ( -1 * (((x - b)**2) / (2*(c**2))) )
    
    return y

def logistic(x, L, k, x0):
    
    f_x = L / (1 + np.exp(-1*k * (x - x0)))
    
    return f_x

def gaussian_logistic_mix(x, a, b, c, L, k, x0):
    
    gaussian_component = gaussian(x, a, b, c)
    logistic_component = logistic(x, L, k, x0)
    y = gaussian_component - logistic_component
    
    return y

def generate_initial_params(ydata, log_mid_guess, log_sup_guess):
    
    p0 = []
    p0.append(np.max(ydata))               # height of the gaussian
    p0.append(len(ydata) / 2)              # center of the gaussian
    p0.append(len(ydata))                  # width of the gaussian
    p0.append(log_sup_guess)               # supremum of the logistic (intensity loss in stripe)
    p0.append(0.05)                        # steepness of the logistic
    p0.append(len(ydata) * log_mid_guess)  # midpoint of the logistic
    
    return p0

def fit_gaussian_logistic_model(ydata, p0, flip_logistic=True, debug_plot=False):
    
    xdata = list(range(len(ydata)))
    
    if flip_logistic:
        p0[4] = -1 * p0[4]
    
    popt, pcov = curve_fit(gaussian_logistic_mix, xdata, ydata, p0=p0, maxfev=int(1e6))
    
    gaussian_fit = gaussian(xdata, popt[0], popt[1], popt[2])
    logistic_fit = logistic(xdata, popt[3], popt[4], popt[5])
    combined_fit = gaussian_logistic_mix(xdata, *popt)
    
    if debug_plot:
        
        fig, axes = plt.subplots(1, 3)
        
        axes[0].plot(xdata, ydata, label='Intensity profile')
        axes[0].plot(xdata, combined_fit, label='Fitted intensity profile')
        axes[1].plot(xdata, ydata, label='Intensity profile')
        axes[1].plot(xdata, gaussian_fit, label='Fitted intensity profile')
        axes[2].plot(xdata, logistic_fit)
        
        axes[0].set_title("Combined model fit")
        axes[1].set_title("Gaussian fit")
        axes[2].set_title("Logistic fit")
        
        axes[0].legend()
        axes[1].legend()
        
        axes[0].set_box_aspect(1)
        axes[1].set_box_aspect(1)
        axes[2].set_box_aspect(1)
        
    return gaussian_fit, logistic_fit, combined_fit

def construct_stripes_flatfield(img, x_logistic_fit, y_logistic_fit):
    
    n_rows, n_cols = img.shape
    stripes_flatfield = np.ones((n_rows, n_cols), dtype=np.float64)
    
    if x_logistic_fit is not None:
        x_logistic_fit = x_logistic_fit.reshape(1, -1)
        stripes_flatfield *= np.tile(1 - x_logistic_fit, (n_rows, 1))
        
    if y_logistic_fit is not None:
        y_logistic_fit = y_logistic_fit.reshape(-1, 1)
        stripes_flatfield *= np.tile(1 - y_logistic_fit, (1, n_cols))
        
    return stripes_flatfield
    
def apply_flatfield_correction(images, flatfield, debug_plot=False):
    
    raw_flatfield = compute_median_flatfield(images) 
    corrected_images = [im / flatfield for im in images]
    corrected_flatfield = compute_median_flatfield(corrected_images)
    
    if debug_plot:
        
        vmin = np.min((np.min(raw_flatfield), np.min(corrected_flatfield)))
        vmax = np.max((np.max(raw_flatfield), np.max(corrected_flatfield)))
        
        fig, axes = plt.subplots(1, 3)
        
        axes[0].imshow(raw_flatfield, vmin=vmin, vmax=vmax)
        axes[1].imshow(flatfield, vmin=vmin, vmax=vmax)
        axes[2].imshow(corrected_flatfield, vmin=vmin, vmax=vmax)
        
        axes[0].set_title("Original flatfield")
        axes[1].set_title("Flatfield correction")
        axes[2].set_title("New flatfield")
    
    return corrected_images, corrected_flatfield  

def build_stripe_flatfield(images, 
                           log_mid_guess,
                           log_sup_guess,
                           correct_x_stripe=True, 
                           correct_y_stripe=True, 
                           flip_logistic=True, 
                           debug_plot=False):
    
    flatfield = compute_median_flatfield(images)
    x_profile, y_profile = compute_intensity_profiles(flatfield, 
                                                      debug_plot=debug_plot)
    
    if correct_x_stripe:
        p0 = generate_initial_params(x_profile, log_mid_guess, log_sup_guess)
        gaussian_fit, x_logistic_fit, combined_fit = fit_gaussian_logistic_model(x_profile, 
                                                                                 p0, 
                                                                                 flip_logistic=flip_logistic, 
                                                                                 debug_plot=debug_plot)
    else:
        x_logistic_fit = None
    
    if correct_y_stripe:
        p0 = generate_initial_params(y_profile, log_mid_guess, log_sup_guess)
        gaussian_fit, y_logistic_fit, combined_fit = fit_gaussian_logistic_model(y_profile, 
                                                                                 p0, 
                                                                                 flip_logistic=flip_logistic, 
                                                                                 debug_plot=debug_plot)
    else:
        y_logistic_fit = None
    
    stripe_flatfield = construct_stripes_flatfield(flatfield, x_logistic_fit, y_logistic_fit)
    # corrected_images, corrected_flatfield = apply_flatfield_correction(images, stripe_flatfield, debug_plot=debug_plot)
    
    return stripe_flatfield

def correct_bleach(corrected_image, position_indices, flatfield, channel_idx, desc=None):
    
    for i in tqdm(position_indices, desc=desc):
        for j in range(n_planes):
            corrected_image[i,j,channel_idx,:,:] = np.uint16(np.round(corrected_image[i,j,channel_idx,:,:] / flatfield))
    
    return corrected_image

def extract_images(image, mask, channel_idx):
    
    im = image[mask,:,channel_idx,:,:]
    t, z, y, x = im.shape
    im = im.reshape(t*z, y, x)
    im = list(im)
    
    return im

def correct_channel(image, channel_idx, log_mid_guess, log_sup_guess, 
                    correct_stripes=True, debug_plot=False):
    
    if correct_stripes:

        print('Fitting left stripe flatfield...')
        left_bleach_im = extract_images(image, left_bleach_mask, channel_idx)
        left_stripe_flatfield = build_stripe_flatfield(left_bleach_im,
                                                       log_mid_guess[0],
                                                       log_sup_guess[0],
                                                       correct_x_stripe=True, 
                                                       correct_y_stripe=False,
                                                       flip_logistic=True,
                                                       debug_plot=debug_plot)
    
        print('Fitting top stripe flatfield...')
        top_bleach_im = extract_images(image, top_bleach_mask, channel_idx)
        top_stripe_flatfield = build_stripe_flatfield(top_bleach_im, 
                                                      log_mid_guess[1],
                                                      log_sup_guess[1],
                                                      correct_x_stripe=False, 
                                                      correct_y_stripe=True,
                                                      flip_logistic=True,
                                                      debug_plot=debug_plot)
    
    
        print('Fitting right stripe flatfield...')
        right_bleach_im = extract_images(image, right_bleach_mask, channel_idx)
        right_stripe_flatfield = build_stripe_flatfield(right_bleach_im,
                                                        log_mid_guess[2],
                                                        log_sup_guess[2],
                                                        correct_x_stripe=True, 
                                                        correct_y_stripe=False,
                                                        flip_logistic=False,
                                                        debug_plot=debug_plot)
        
        flatfield_2 = left_stripe_flatfield
        flatfield_3 = top_stripe_flatfield
        flatfield_4 = top_stripe_flatfield * right_stripe_flatfield
        flatfield_5 = top_stripe_flatfield * left_stripe_flatfield
        
        image = correct_bleach(image, indices_2, flatfield_2, channel_idx, desc='Correcting flatfield 2 images...')
        image = correct_bleach(image, indices_3, flatfield_3, channel_idx, desc='Correcting flatfield 3 images...')
        image = correct_bleach(image, indices_4, flatfield_4, channel_idx, desc='Correcting flatfield 4 images...')
        image = correct_bleach(image, indices_5, flatfield_5, channel_idx, desc='Correcting flatfield 5 images...')

    print('Fitting BaSiC object...')
    basic_obj = fit_channel_basic(image, channel_idx)
    image = correct_channel_basic(image, basic_obj, channel_idx)
    
    return image

def correct_image(image):
    
    print('Making copy of image...')
    corrected_image = image.copy()
    
    correct_stripes = [True, False, False, False, False]
    log_mid_guess = [0.17, 0.17, 0.83]
    log_sup_guess = [0.07, 0.07, 0.07]
    
    for i in range(n_chans):
        print(f'Correcting image channel index {i}...')
        corrected_image = correct_channel(corrected_image, i, log_mid_guess, 
                                          log_sup_guess, correct_stripes=correct_stripes[i])
        
    return corrected_image

#%%

corrected_image = correct_image(my_array)

#%% export to OME-TIFF
# nightmarish code to export the shading corrected file as an OME-TIFF with all
# the requisite metadata inherited from the original nd2 file 

from ome_types.model import OME, Image, Pixels, Plane, Channel, Pixels_DimensionOrder, TiffData

print('Saving image file...')

sizes = myfile.sizes
voxel = myfile.voxel_size()

positions = []
for loop in myfile.experiment:
    if loop.type == "XYPosLoop":
        positions = [(p.stagePositionUm.x, p.stagePositionUm.y, p.stagePositionUm.z) for p in loop.parameters.points]

assert len(positions) == sizes["P"]
p_axis = list(sizes).index("P")
n_c = sizes.get("C", 1)
n_z = sizes.get("Z", 1)

images = []
first_ifd = 0
planes_per_series = n_c * n_z

for i, (x, y, z) in enumerate(positions):
 
    id=f"Image:{i}",
    name=f"tile_{i}",
    pixels=Pixels(
        id=f"Pixels:{i}",
        dimension_order=Pixels_DimensionOrder.XYCZT,
        size_x=sizes["X"], size_y=sizes["Y"],
        size_c=n_c, size_z=n_z, size_t=1,
        type="uint16",
        physical_size_x=voxel.x, physical_size_y=voxel.y,
        physical_size_z=voxel.z if n_z > 1 else None,
        channels=[Channel(id=f"Channel:{i}:{c}", samples_per_pixel=1) for c in range(n_c)],
        planes=[
            Plane(the_c=c, the_z=zi, the_t=0, position_x=x, position_y=y,
                  position_z=z + zi * voxel.z if voxel.z else z)
            for zi in range(n_z)
            for c in range(n_c)
        ],
        tiff_data_blocks=[
            TiffData(ifd=first_ifd, plane_count=planes_per_series)
        ],
    )
    
    images.append(Image(id=f"Image:{i}", name=f"tile_{i}", pixels=pixels))
    first_ifd += planes_per_series   # advance cumulative IFD offset for next series
 
ome = OME(images=images)

tiles = [np.take(corrected_image, i, axis=p_axis) for i in range(sizes["P"])]

output_path = os.path.join(parent, stem + ".tiff")

with tifffile.TiffWriter(output_path, bigtiff=True) as tif:
    for i, tile in enumerate(tiles):
        tif.write(
            tile,
            photometric="minisblack",
            compression="zlib",
            contiguous=False,           # force a new series, don't chain onto previous
            metadata={"axes": "ZCYX"},  # explicit axes, prevents shape-based auto-merge
            description=ome.to_xml() if i == 0 else None,  # OME-XML goes on IFD 0 only
        )

#%%

def plot_flatfields():
    
    fig, axes = plt.subplots(2, 2)
    
    vmin = np.min((np.min(flatfield_2),
                  np.min(flatfield_3),
                  np.min(flatfield_4),
                  np.min(flatfield_5)))
    
    vmax = np.max((np.max(flatfield_2),
                  np.max(flatfield_3),
                  np.max(flatfield_4),
                  np.max(flatfield_5)))

    im2 = axes[0,0].imshow(flatfield_2, vmin=vmin, vmax=vmax)
    im3 = axes[0,1].imshow(flatfield_3, vmin=vmin, vmax=vmax)
    im4 = axes[1,0].imshow(flatfield_4, vmin=vmin, vmax=vmax)
    im5 = axes[1,1].imshow(flatfield_5, vmin=vmin, vmax=vmax)
    
    axes[0,0].set_title('Flatfield 2')
    axes[0,1].set_title('Flatfield 3')
    axes[1,0].set_title('Flatfield 4')
    axes[1,1].set_title('Flatfield 5')
    
    cbar2 = fig.colorbar(im2, ax=axes[0,0])
    cbar2.set_label('Relative intensity', rotation=270, labelpad=15, fontsize=11)
    cbar3 = fig.colorbar(im3, ax=axes[0,1])
    cbar3.set_label('Relative intensity', rotation=270, labelpad=15, fontsize=11)
    cbar4 = fig.colorbar(im4, ax=axes[1,0])
    cbar4.set_label('Relative intensity', rotation=270, labelpad=15, fontsize=11)
    cbar5 = fig.colorbar(im5, ax=axes[1,1])
    cbar5.set_label('Relative intensity', rotation=270, labelpad=15, fontsize=11)
    
if debug_plot:
    plot_flatfields()



    

        


