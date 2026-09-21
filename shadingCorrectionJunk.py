# -*- coding: utf-8 -*-
"""
Created on Tue Aug 11 19:22:17 2026

@author: jpv88
"""

#%%

left_corner_bleach_im = extract_images(corrected_image, left_corner_bleach_mask, 0)
left_corner_flatfield = build_stripe_flatfield(left_corner_bleach_im, 
                                   correct_x_stripe=False, 
                                   correct_y_stripe=True,
                                   flip_logistic=True,
                                   x_percentile=10,
                                   y_percentile=10,
                                   debug_plot=True)

lo, hi = left_stripe_flatfield.min(), left_stripe_flatfield.max()
left_stripe_flatfield_inverted = (hi - left_stripe_flatfield) / (hi - lo)

lo, hi = top_stripe_flatfield.min(), top_stripe_flatfield.max()
top_stripe_flatfield_inverted = (hi - top_stripe_flatfield) / (hi - lo)

lo, hi = right_stripe_flatfield.min(), right_stripe_flatfield.max()
right_stripe_flatfield_inverted = (hi - right_stripe_flatfield) / (hi - lo)

corner_mask = (top_stripe_flatfield_inverted * left_stripe_flatfield_inverted)
corner_mask = corner_mask * left_corner_flatfield.min()
corner_mask = 1 - corner_mask * (1 - left_corner_flatfield.min()) / left_corner_flatfield.min()
left_corner_flatfield = corner_mask

right_corner_bleach_im = extract_images(corrected_image, right_corner_bleach_mask, 0)
right_corner_flatfield = build_stripe_flatfield(right_corner_bleach_im, 
                                   correct_x_stripe=False, 
                                   correct_y_stripe=True,
                                   flip_logistic=True,
                                   x_percentile=90,
                                   y_percentile=90,
                                   debug_plot=True)

corner_mask = (top_stripe_flatfield_inverted * right_stripe_flatfield_inverted)
corner_mask = corner_mask * right_corner_flatfield.min()
corner_mask = 1 - corner_mask * (1 - right_corner_flatfield.min()) / right_corner_flatfield.min()
right_corner_flatfield = corner_mask

#%%

flatfield_2 = left_stripe_flatfield
flatfield_3 = top_stripe_flatfield * left_corner_flatfield
flatfield_4 = top_stripe_flatfield * right_stripe_flatfield * left_corner_flatfield * right_corner_flatfield
flatfield_5 = top_stripe_flatfield * right_corner_flatfield * right_stripe_flatfield
flatfield_6 = top_stripe_flatfield * right_corner_flatfield
flatfield_7 = top_stripe_flatfield * left_stripe_flatfield * left_corner_flatfield * right_corner_flatfield
flatfield_8 = top_stripe_flatfield * left_stripe_flatfield * left_corner_flatfield

def plot_flatfields():
    
    fig, axes = plt.subplots(2, 4)
    
    flatfield_1 = np.ones(flatfield_2.shape)
    
    vmin = np.min((np.min(flatfield_1),
                  np.min(flatfield_2),
                  np.min(flatfield_3),
                  np.min(flatfield_4),
                  np.min(flatfield_5),
                  np.min(flatfield_6),
                  np.min(flatfield_7),
                  np.min(flatfield_8)))
    
    vmax = np.max((np.max(flatfield_1),
                  np.max(flatfield_2),
                  np.max(flatfield_3),
                  np.max(flatfield_4),
                  np.max(flatfield_5),
                  np.max(flatfield_6),
                  np.max(flatfield_7),
                  np.max(flatfield_8)))
    
    axes[0,0].imshow(flatfield_1, vmin=vmin, vmax=vmax)
    axes[0,1].imshow(flatfield_2, vmin=vmin, vmax=vmax)
    axes[0,2].imshow(flatfield_3, vmin=vmin, vmax=vmax)
    axes[0,3].imshow(flatfield_4, vmin=vmin, vmax=vmax)
    axes[1,0].imshow(flatfield_5, vmin=vmin, vmax=vmax)
    axes[1,1].imshow(flatfield_6, vmin=vmin, vmax=vmax)
    axes[1,2].imshow(flatfield_7, vmin=vmin, vmax=vmax)
    axes[1,3].imshow(flatfield_8, vmin=vmin, vmax=vmax)
    
plot_flatfields()
    

#%%

corrected_image = my_array.copy()

channel_idx = 0

corrected_image = correct_bleach(corrected_image, indices_2, flatfield_2, channel_idx)
corrected_image = correct_bleach(corrected_image, indices_3, flatfield_3, channel_idx)
corrected_image = correct_bleach(corrected_image, indices_4, flatfield_4, channel_idx)
corrected_image = correct_bleach(corrected_image, indices_5, flatfield_5, channel_idx)
corrected_image = correct_bleach(corrected_image, indices_6, flatfield_6, channel_idx)
corrected_image = correct_bleach(corrected_image, indices_7, flatfield_7, channel_idx)
corrected_image = correct_bleach(corrected_image, indices_8, flatfield_8, channel_idx)

def correct_left_bleach(corrected_image, flatfield, channel_idx):
    
    for i in indices_left_bleach:
        for j in range(n_planes):
            corrected_image[i,j,channel_idx,:,:] = np.uint16(np.round(corrected_image[i,j,channel_idx,:,:] / flatfield))
    
    return corrected_image

def correct_top_bleach(corrected_image, flatfield, channel_idx):
    
    for i in indices_top_bleach:
        for j in range(n_planes):
            corrected_image[i,j,channel_idx,:,:] = np.uint16(np.round(corrected_image[i,j,channel_idx,:,:] / flatfield))
            
        return corrected_image

def correct_right_bleach(corrected_image, flatfield, channel_idx):
    
    for i in indices_right_bleach:
        for j in range(n_planes):
            corrected_image[i,j,channel_idx,:,:] = np.uint16(np.round(corrected_image[i,j,channel_idx,:,:] / flatfield))
    
    return corrected_image

#%%

def plot_flatfields(basic_objs):
    
    fig = plt.figure(figsize=(10, 6))
    gs = GridSpec(2, 6, figure=fig)  # use 6 columns for easy divisibility
    
    ax1 = fig.add_subplot(gs[0, 0:2])
    ax2 = fig.add_subplot(gs[0, 2:4])
    ax3 = fig.add_subplot(gs[0, 4:6])
    ax4 = fig.add_subplot(gs[1, 1:3])
    ax5 = fig.add_subplot(gs[1, 3:5])

    im1 = ax1.imshow(basic_objs[0].flatfield)
    fig.colorbar(im1, ax=ax1)
    ax1.set_title('Photobleaching type 1 (no bleach)')
    
    im2 = ax2.imshow(basic_objs[1].flatfield)
    fig.colorbar(im2, ax=ax2)
    ax2.set_title('Photobleaching type 2 (left bleach only)')
    
    im3 = ax3.imshow(basic_objs[2].flatfield)
    fig.colorbar(im3, ax=ax3)
    ax3.set_title('Photobleaching type 3 (top bleach only)')
    
    im4 = ax4.imshow(basic_objs[3].flatfield)
    fig.colorbar(im4, ax=ax4)
    ax4.set_title('Photobleaching type 4 (top bleach and right bleach)')
    
    im5 = ax5.imshow(basic_objs[4].flatfield)
    fig.colorbar(im5, ax=ax5)
    ax5.set_title('Photobleaching type 5 (top bleach and left bleach)')
    
    plt.tight_layout()
    plt.show()
    
plot_flatfields(basic_objs)

#%% 

# photobleaching types
# type 1: no bleach
# type 2: left bleach only
# type 3: top bleach and top left corner bleach
# type 4: top bleach, top left corner bleach, top right corner bleach, right bleach
# type 5: top bleach, top right corner bleach, right bleach
# type 6: top bleach and top right corner bleach
# type 7: top bleach, top left corner bleach, top right corner bleach, left bleach
# type 8: top bleach, top left corner bleach, left bleach

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
        
    # even row, first column
    elif rows[i] % 2 == 0 and cols[i] == 1:
        bleach_types.append(5)
    
    # even row, not first or last column
    elif rows[i] % 2 == 0:
        bleach_types.append(4)
    
    # odd row, first column
    elif rows[i] % 2 != 0 and cols[i] == 1:
        bleach_types.append(6)
        
    # odd row, last column
    elif rows[i] % 2 != 0 and cols[i] == n_cols:
        bleach_types.append(8)
    
    # odd row, not first or last column
    elif rows[i] % 2 != 0:
        bleach_types.append(7)
        
mask_1 = (np.array(bleach_types) == 1)
mask_2 = (np.array(bleach_types) == 2)
mask_3 = (np.array(bleach_types) == 3)
mask_4 = (np.array(bleach_types) == 4)
mask_5 = (np.array(bleach_types) == 5)
mask_6 = (np.array(bleach_types) == 6)
mask_7 = (np.array(bleach_types) == 7)
mask_8 = (np.array(bleach_types) == 8)

indices_1 = np.arange(n_pos)[mask_1]
indices_2 = np.arange(n_pos)[mask_2]
indices_3 = np.arange(n_pos)[mask_3]
indices_4 = np.arange(n_pos)[mask_4]
indices_5 = np.arange(n_pos)[mask_5]
indices_6 = np.arange(n_pos)[mask_6]
indices_7 = np.arange(n_pos)[mask_7]
indices_8 = np.arange(n_pos)[mask_8]

left_bleach_mask = mask_2 | mask_7 | mask_8
indices_left_bleach = np.arange(n_pos)[left_bleach_mask]

right_bleach_mask = mask_4 | mask_5
indices_right_bleach = np.arange(n_pos)[right_bleach_mask]

top_bleach_mask = mask_3 | mask_4 | mask_5 | mask_6 | mask_7 | mask_8
indices_top_bleach = np.arange(n_pos)[top_bleach_mask]

left_corner_bleach_mask = mask_3 | mask_4 | mask_7 | mask_8
indices_left_corner_bleach = np.arange(n_pos)[left_corner_bleach_mask]

right_corner_bleach_mask = mask_4 | mask_5 | mask_6 | mask_7
indices_right_corner_bleach = np.arange(n_pos)[right_corner_bleach_mask]

#%%

# # fit basic object for a single bleaching type and channel
# def fit_basic_bleach(basic_obj, images):
    
#     t, z, y, x = images.shape
#     images = images.reshape(t*z, y, x)
    
#     basic_obj.fit(images)
    
#     return basic_obj

# # fit basic object for all bleaching types in a channel
# def fit_basic_channel(image, channel_idx):
    
#     images_1 = image[mask_1,:,channel_idx,:,:]
#     images_2 = image[mask_2,:,channel_idx,:,:]
#     images_3 = image[mask_3,:,channel_idx,:,:]
#     images_4 = image[mask_4,:,channel_idx,:,:]
#     images_5 = image[mask_5,:,channel_idx,:,:]
    
#     basic_1 = BaSiC(get_darkfield=False, fitting_mode='ladmap', smoothness_flatfield=smoothness_flatfield, sort_intensity=True)
#     basic_2 = BaSiC(get_darkfield=False, fitting_mode='ladmap', smoothness_flatfield=smoothness_flatfield, sort_intensity=True)
#     basic_3 = BaSiC(get_darkfield=False, fitting_mode='ladmap', smoothness_flatfield=smoothness_flatfield, sort_intensity=True)
#     basic_4 = BaSiC(get_darkfield=False, fitting_mode='ladmap', smoothness_flatfield=smoothness_flatfield, sort_intensity=True)
#     basic_5 = BaSiC(get_darkfield=False, fitting_mode='ladmap', smoothness_flatfield=smoothness_flatfield, sort_intensity=True)

#     with tqdm(total=5, desc=f"Calculating shading model for bleaching type 1, channel index {channel_idx}") as pbar:
        
#         basic_1 = fit_basic_bleach(basic_1, images_1)
#         pbar.update(1)
        
#         pbar.set_description(f"Calculating shading model for bleaching type 2, channel index {channel_idx}")
#         basic_2 = fit_basic_bleach(basic_2, images_2)
#         pbar.update(1)
        
#         pbar.set_description(f"Calculating shading model for bleaching type 3, channel index {channel_idx}")
#         basic_3 = fit_basic_bleach(basic_3, images_3)
#         pbar.update(1)
        
#         pbar.set_description(f"Calculating shading model for bleaching type 4, channel index {channel_idx}")
#         basic_4 = fit_basic_bleach(basic_4, images_4)
#         pbar.update(1)
        
#         pbar.set_description(f"Calculating shading model for bleaching type 5, channel index {channel_idx}")
#         basic_5 = fit_basic_bleach(basic_5, images_5)
#         pbar.update(1)
    
#     return [basic_1, basic_2, basic_3, basic_4, basic_5]

# def correct_channel(image, corrected_image, basic_objs, channel_idx):
    
#     with tqdm(total=5, desc=f"Correcting images for bleaching type 1, channel index {channel_idx}") as pbar:
    
#         for i in indices_1:
#             for j in range(n_planes):
#                 input_image = image[i,j,channel_idx,:,:]
#                 input_image = np.expand_dims(input_image, axis=0)
#                 corrected_image[i,j,channel_idx,:,:] = basic_objs[0].transform(input_image, use_tqdm=False)
#         pbar.update(1)
        
#         pbar.set_description(f"Correcting images for bleaching type 2, channel index {channel_idx}")
#         for i in indices_2:
#             for j in range(n_planes):
#                 input_image = image[i,j,channel_idx,:,:]
#                 input_image = np.expand_dims(input_image, axis=0)
#                 corrected_image[i,j,channel_idx,:,:] = basic_objs[1].transform(input_image, use_tqdm=False)
#         pbar.update(1)
        
#         pbar.set_description(f"Correcting images for bleaching type 3, channel index {channel_idx}")
#         for i in indices_3:
#             for j in range(n_planes):
#                 input_image = image[i,j,channel_idx,:,:]
#                 input_image = np.expand_dims(input_image, axis=0)
#                 corrected_image[i,j,channel_idx,:,:] = basic_objs[2].transform(input_image, use_tqdm=False)
#         pbar.update(1)
        
#         pbar.set_description(f"Correcting images for bleaching type 4, channel index {channel_idx}")
#         for i in indices_4:
#             for j in range(n_planes):
#                 input_image = image[i,j,channel_idx,:,:]
#                 input_image = np.expand_dims(input_image, axis=0)
#                 corrected_image[i,j,channel_idx,:,:] = basic_objs[3].transform(input_image, use_tqdm=False)
#         pbar.update(1)
              
#         pbar.set_description(f"Correcting images for bleaching type 5, channel index {channel_idx}")
#         for i in indices_5:
#             for j in range(n_planes):
#                 input_image = image[i,j,channel_idx,:,:]
#                 input_image = np.expand_dims(input_image, axis=0)
#                 corrected_image[i,j,channel_idx,:,:] = basic_objs[4].transform(input_image, use_tqdm=False)
#         pbar.update(1)
    
#     return corrected_image

#%%

channel_idx = 0
im = extract_images(my_array, left_bleach_mask, channel_idx)
flatfield = compute_median_flatfield(im)
plt.imshow(flatfield)


#%%

mask = np.ones(n_pos, dtype=bool)

im = extract_images(corrected_image, mask, 2)
flatfield = compute_median_flatfield(im)


#%%

# corrected_image = np.zeros(myfile.shape, dtype=np.uint16)

# for i in range(n_chans):
#     basic_objs = fit_basic_channel(my_array, i)
#     corrected_image = correct_channel(my_array, corrected_image, basic_objs, i)
    

corrected_image = np.zeros(myfile.shape, dtype=np.uint16)

for i in range(n_chans):
    basic_obj = fit_basic_channel(my_array, i)
    corrected_image = correct_channel(my_array, corrected_image, basic_obj, i)
    
#%%


from scipy.ndimage import gaussian_filter

test_chan = my_array[mask_2,:,0,:,:]
t, z, y, x = test_chan.shape
test_chan = test_chan.reshape(t*z, y, x)

test_chan = list(test_chan)

#%%

flatfield = compute_median_flatfield(test_chan)
x_profile, y_profile = compute_intensity_profiles(flatfield, debug_plot=True)
p0 = generate_initial_params(x_profile)
gaussian_fit, logistic_fit, combined_fit = fit_gaussian_logistic_model(x_profile, p0, debug_plot=True)
stripes_flatfield = construct_stripes_flatfield(flatfield, logistic_fit, None)
corrected_images, corrected_flatfield = apply_flatfield_correction(test_chan, stripes_flatfield, debug_plot=True)

#%%

row_profile = np.median(typical, axis=1)
row_profile_smooth = uniform_filter1d(row_profile, size=10)
row_profile_smooth = row_profile_smooth.reshape(-1, 1)

col_profile = np.median(typical, axis=0)
col_profile_smooth = uniform_filter1d(col_profile, size=10)

# test_chan_corr = [a / row_profile_smooth for a in test_chan]
# test_chan_corr = [a / col_profile_smooth for a in test_chan_corr]

# # 1. normalize each image by its own brightness
# normed = [img / np.median(img) for img in test_chan_corr]

# # 2. median-stack to get the "typical" fixed-pattern image
# stack = np.stack(normed, axis=0)
# typical = np.median(stack, axis=0)

#%%

basic_obj = fit_basic_list(test_chan)
plt.imshow(basic_obj.flatfield)

#%%



def f_gauss(x, a, b, c):
    
    f_x = a * np.exp ( -1 * ( (x - b)**2 / (2 * c**2) ) )
    
    return f_x

def f_reverse_sigmoid(x, L, k, midpoint):
    
    f_x = L / (1 + np.exp(k*(x - midpoint)))
    
    return f_x

def f(x, a, b, c, L, k, midpoint):
    
    y_gauss = f_gauss(x, a, b, c)
    y_reverse_sigmoid = f_reverse_sigmoid(x, L, k, midpoint)
    f_x = y_gauss - y_reverse_sigmoid
    
    return f_x

target_vec = row_profile_smooth
xdata = list(range(len(target_vec)))
ydata = np.squeeze(target_vec)

p0 = []
p0.append(np.max(target_vec))
p0.append(len(target_vec)/2)
p0.append(1000)
p0.append(0.05)
p0.append(0.05)
p0.append(300)
popt, pcov = curve_fit(f, xdata, ydata, p0=p0)
fit = f(xdata, *popt)
fit_reverse_sigmoid = f_reverse_sigmoid(xdata, popt[3], popt[4], popt[5])
fit_reverse_sigmoid = fit_reverse_sigmoid.reshape(-1, 1)

fig, ax = plt.subplots()
plt.plot(xdata, ydata)
plt.plot(xdata, fit)

bleach_flatfield = np.tile(1 - fit_reverse_sigmoid, (1, 976))

#%%

test_chan_corr = [a / bleach_flatfield for a in test_chan]

# # 1. normalize each image by its own brightness
normed = [img / np.median(img) for img in test_chan_corr]

# # 2. median-stack to get the "typical" fixed-pattern image
stack = np.stack(normed, axis=0)
typical = np.median(stack, axis=0)

# fit basic object for a list of images
def fit_basic_list(image):
    
    image = np.stack(image, axis=0)
    
    basic_obj = BaSiC(get_darkfield=False, fitting_mode='ladmap', sort_intensity=True)
    basic_obj.fit(image)
        
    return basic_obj