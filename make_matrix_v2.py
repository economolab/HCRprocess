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

data_dir = 'D:\\2026-09-12_FN_SNr_2\\post\\core_output\\s01L'

uniq_id, post_dir, exp_name = matrix_utils.derive_params(data_dir)
data_f, reg_f, bin_f, mask_f, mask_im_f = matrix_utils.build_paths(data_dir, 
                                                                   uniq_id, 
                                                                   post_dir, 
                                                                   exp_name,
                                                                   mode='curated')

remove_b = False

#%% load everything

mask_df = pd.read_csv(os.path.join(data_dir,mask_f))
bin_df = pd.read_csv(os.path.join(data_dir,bin_f))
cell_df = matrix_utils.build_cell_df(mask_df, bin_df)

data_f = matrix_utils.remove_label_channels(data_f)
quant_genes = [f[3:-7] for f in data_f]

if remove_b:
    cell_df = cell_df.loc[:, ~cell_df.columns.str.contains('_b')]
else:
    quant_genes.remove('slc17a6')
    quant_genes.remove('slc32a1')
    quant_genes.remove('rph3a')
    quant_genes.remove('slc5a7')
    data_f = data_f[data_f != 'r1_slc17a6488.tif']
    data_f = data_f[data_f != 'r2_slc32a1488.tif']
    data_f = data_f[data_f != 'r3_rph3a488.tif']
    data_f = data_f[data_f != 'r4_slc5a7488.tif']
    
for quant_gene in quant_genes:
    cell_df.insert(cell_df.shape[1],quant_gene,np.zeros(cell_df.shape[0]))
    
cell_df['Cell ID'] = cell_df['Cell ID'].astype(str)
for i in range(cell_df.shape[0]):
    cell_df.loc[i,"Cell ID"] = exp_name + '_' + uniq_id + '_' + str(i+1)
    
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

    
#%% load images

quant_gene_ims = []

for i in tqdm(range(len(data_f)),desc='Loading images...'):
    im = tifffile.imread(os.path.join(data_dir,data_f[i]))
    im = np.transpose(im,[1,2,0])
    im[np.isnan(im)] = 0
    quant_gene_ims.append(im)
    
quant_gene_ims = dict(zip(quant_genes, quant_gene_ims))
    
masks_im = tifffile.imread(mask_im_f)
masks_im = np.transpose(masks_im,[1,2,0])

pp_masks_im = filter_data.filt_masks_im_pp(masks_im, cell_df)
row_slice, col_slice = img_utils.make_crop_slices(pp_masks_im)

#%%

from bigfish.stack import log_filter
from bigfish.detection.spot_detection import local_maximum_detection, spots_thresholding
from kneed import KneeLocator

from scipy.interpolate import UnivariateSpline

def set_thresh(image, fluo, inlier_planes, min_spots=int(100), max_spots=int(1e5), 
               num_evals=100, gene=None, save_threshs=True, plot_fit=False):
    
    psf_sigma = img_utils.calc_psf_sigma(fluo)
    voxel_size = 325
    sigma = psf_sigma / voxel_size
    minimum_distance = sigma
     
    im_filtered = np.zeros_like(image)
    n_planes = image.shape[2]
     
    for k in tqdm(range(n_planes)):
        im_filtered[:,:,k] = log_filter(image[:,:,k], sigma)
    
    first_thresh = []
    first_spots = []
    for k in tqdm(range(n_planes)):
        
        spots_mask = local_maximum_detection(im_filtered[:,:,k], minimum_distance)
        min_px = np.min(im_filtered[:,:,k])
        max_px = np.max(im_filtered[:,:,k])
        threshs = np.linspace(min_px, max_px, num_evals)
        
        num_spots = []
        for thresh in tqdm(threshs):
            spots, _ = spots_thresholding(im_filtered[:,:,k], spots_mask, thresh, remove_duplicate=False)
            num_spots.append(len(spots))
            
        thresh_spots = np.column_stack((threshs, num_spots))
        _, idx = np.unique(thresh_spots[::-1, 1], return_index=True)
        orig_idx = thresh_spots.shape[0] - 1 - idx
        thresh_spots = thresh_spots[np.sort(orig_idx), :]
        
        thresh_spots = thresh_spots[thresh_spots[:,1] >= min_spots, :]
        thresh_spots = thresh_spots[thresh_spots[:,1] <= max_spots, :]
        
        kl = KneeLocator(thresh_spots[:,0], thresh_spots[:,1], curve="convex", direction="decreasing")
        first_thresh.append(kl.knee)
        
        spots, _ = spots_thresholding(im_filtered[:,:,k], spots_mask, kl.knee, remove_duplicate=False)
        first_spots.append(len(spots))
        
    # spots_inlier = np.array(final_spots)[inlier_planes]
    threshs_inlier = np.array(first_thresh)[inlier_planes]
    
    z_plane = np.arange(n_planes)
    z_plane_inlier = z_plane[inlier_planes]
    smooth_thresh = matrix_utils.fit_threshold_curve(z_plane_inlier, np.array(threshs_inlier),
                                                     plot_fit=plot_fit, predict_z_plane=z_plane)
    
    # spots_std = np.std(spots_inlier)
    # spots_mean = np.mean(spots_inlier)
    
    # spots_range = [spots_mean - 3*spots_std, spots_mean + 3*spots_std]
    # spots_range = np.clip(spots_range, a_min=min_spots, a_max=max_spots)
    
    smooth_thresh_ptp = np.ptp(smooth_thresh)
    threshs_range = [min(smooth_thresh) - (smooth_thresh_ptp/2), 
                     max(smooth_thresh) + (smooth_thresh_ptp/2)]
    
    final_thresh = []
    final_spots = []
    for k in tqdm(range(n_planes)):
        
        spots_mask = local_maximum_detection(im_filtered[:,:,k], minimum_distance)
        min_px = np.min(im_filtered[:,:,k])
        max_px = np.max(im_filtered[:,:,k])
        # threshs_range_plane = [smooth_thresh[k] - (smooth_thresh_ptp/2), smooth_thresh[k] + (smooth_thresh_ptp/2)]
        threshs_range = np.clip(threshs_range, a_min=min_px, a_max=max_px)
        threshs = np.linspace(threshs_range[0], threshs_range[1], num_evals)
        
        num_spots = []
        for thresh in tqdm(threshs):
            spots, _ = spots_thresholding(im_filtered[:,:,k], spots_mask, thresh, remove_duplicate=False)
            num_spots.append(len(spots))
            
        thresh_spots = np.column_stack((threshs, num_spots))
        _, idx = np.unique(thresh_spots[::-1, 1], return_index=True)
        orig_idx = thresh_spots.shape[0] - 1 - idx
        thresh_spots = thresh_spots[np.sort(orig_idx), :]
        
        # thresh_spots = thresh_spots[thresh_spots[:,1] >= spots_range[0], :]
        # thresh_spots = thresh_spots[thresh_spots[:,1] <= spots_range[1], :]
        
        thresh_spots = thresh_spots[thresh_spots[:,1] >= min_spots, :]
        thresh_spots = thresh_spots[thresh_spots[:,1] <= max_spots, :]
        
        kl = KneeLocator(thresh_spots[:,0], thresh_spots[:,1], curve="convex", direction="decreasing")
        final_thresh.append(kl.knee)
        
        spots, _ = spots_thresholding(im_filtered[:,:,k], spots_mask, kl.knee, remove_duplicate=False)
        final_spots.append(len(spots))
        
    z_plane = np.arange(n_planes)
    z_plane_inlier = z_plane[inlier_planes]
    threshs_inlier = np.array(final_thresh)[inlier_planes]
    smooth_thresh = matrix_utils.fit_threshold_curve(z_plane_inlier, np.array(threshs_inlier), 
                                                     plot_fit=plot_fit, predict_z_plane=z_plane)
    
    if save_threshs:
        save_path = os.path.join(data_dir, gene + '_threshs.csv')
        np.savetxt(save_path, smooth_thresh, delimiter=',')
    
    return smooth_thresh
    
    
def count_spots_per_cell(cell_df, gene, fluo, image, pp_masks_im, threshs, return_num_spots=False):
    
    n_planes = image.shape[2]
    
    num_spots = []
    cell_spots_dict = {}
    for cell_id in cell_df['Mask ID'].values:
        cell_spots_dict[cell_id] = 0
        
    for k in tqdm(range(n_planes)):
        
        spots = matrix_utils.count_spots(image[:,:,k], fluo, thresh=threshs[k])
        num_spots.append(len(spots))
        
        masks_at_spots = pp_masks_im[:,:,k][spots[:,0], spots[:,1]]
        unique_masks, unique_counts = np.unique(masks_at_spots, return_counts=True)
        keep = (unique_masks != 0)
        unique_masks, unique_counts = unique_masks[keep], unique_counts[keep]
        
        if len(unique_masks) != 0:
            for idx, mask in enumerate(unique_masks):
                cell_spots_dict[mask] += unique_counts[idx]
        
    cell_df[gene] = list(cell_spots_dict.values())
    
    if return_num_spots:
        return cell_df, num_spots
    else:
        return cell_df
    

def quantify_gene(cell_df, gene, image, pp_masks_im, load_thresh=True, 
                  crop_im=True, row_slice=None, col_slice=None, return_num_spots=False):
    
    fluo = gene_to_fluo[gene]
    
    if crop_im:
        image = image[row_slice, col_slice, :]
        pp_masks_im = pp_masks_im[row_slice, col_slice, :]
    
    image, _ = img_utils.correct_depth_attenuation(image, 
                                                   scale_to_log_filter=True,
                                                   plot_fit=False,
                                                   fluo=fluo)
    
    if load_thresh:
        threshs_path = os.path.join(data_dir, gene + '_threshs.csv')
        if os.path.exists(threshs_path):
            threshs = np.loadtxt(threshs_path, delimiter=',')
        else:
            threshs = set_thresh(image, fluo, gene=gene)
    else:
        threshs = set_thresh(image, fluo, gene=gene)
    
    if return_num_spots:
        cell_df, num_spots = count_spots_per_cell(cell_df, gene, fluo, image, pp_masks_im, threshs, return_num_spots=True)
        return cell_df, num_spots
    else:
        cell_df = count_spots_per_cell(cell_df, gene, fluo, image, pp_masks_im, threshs)
        return cell_df

def calc_noise_floor(im, c=6):
    
    med = np.median(im)
    mad = 1.4826 * np.median(np.abs(im - med))
    floor_im = med + c * mad   # c around 4-6, tune on a deep plane
    
    return floor_im
    

#%%

gene = 'meis2'
crop_im = True

fluo = gene_to_fluo[gene]
image = quant_gene_ims[gene].copy()

if crop_im:
    image = image[row_slice, col_slice, :]

image, inlier_mask = img_utils.correct_depth_attenuation(image, 
                                            scale_to_log_filter=True,
                                            plot_fit=True,
                                            fluo=fluo)

set_thresh(image, gene_to_fluo[gene], inlier_mask, gene=gene, num_evals=100, plot_fit=True)

cell_df, num_spots = quantify_gene(cell_df,
                                   gene,
                                   quant_gene_ims[gene],
                                   pp_masks_im,
                                   crop_im=True,
                                   row_slice=row_slice,
                                   col_slice=col_slice,
                                   return_num_spots=True)

#%% 

# genes_to_quant = ['phox2b', 'ralyl', 'zfhx3']

for gene in quant_genes:
    
    print(gene)
    fluo = gene_to_fluo[gene]
    image = quant_gene_ims[gene].copy()
    image = image[row_slice, col_slice, :]

    image, inlier_mask = img_utils.correct_depth_attenuation(image, 
                                                scale_to_log_filter=True,
                                                plot_fit=False,
                                                fluo=fluo)
    
    set_thresh(image, fluo, inlier_mask, gene=gene, num_evals=100, plot_fit=True)

#%% 

# genes_to_quant = ['phox2b', 'ralyl', 'zfhx3']

num_spots_dict = {}

for gene in quant_genes:
    
    cell_df, num_spots_dict[gene] = quantify_gene(cell_df, gene, quant_gene_ims[gene], pp_masks_im,
                            crop_im=True, row_slice=row_slice, col_slice=col_slice, return_num_spots=True)
    
#%%

import matplotlib.pyplot as plt

for gene in quant_genes:
    
    plt.subplots()
    plt.plot(num_spots_dict[gene])
    plt.title(gene)
    
#%% 1. per-cell, per-plane table

import numpy as np, pandas as pd
import statsmodels.formula.api as smf
import statsmodels.api as sm
import matplotlib.pyplot as plt

gene = 'meis2'
fluo = gene_to_fluo[gene]
im = quant_gene_ims[gene][row_slice, col_slice, :]
masks = pp_masks_im[row_slice, col_slice, :]
im, _ = img_utils.correct_depth_attenuation(im, scale_to_log_filter=True,
                                            plot_fit=False, fluo=fluo)
threshs = np.loadtxt(os.path.join(data_dir, gene + '_threshs.csv'), delimiter=',')

def per_plane_table(im, masks, fluo, threshs, scale=1.0):
    rows = []
    for k in tqdm(range(im.shape[2])):
        spots = matrix_utils.count_spots(im[:,:,k], fluo, thresh=threshs[k]*scale)
        m = masks[:,:,k]
        ids, areas = np.unique(m[m > 0], return_counts=True)
        if len(spots):
            hit = m[spots[:,0], spots[:,1]]
            cnt = pd.Series(hit[hit > 0]).value_counts()
        else:
            cnt = pd.Series(dtype=int)
        for cid, a in zip(ids, areas):
            rows.append((cid, k, a, int(cnt.get(cid, 0))))
    return pd.DataFrame(rows, columns=['mask_id', 'plane', 'area', 'n_spots'])

df = per_plane_table(im, masks, fluo, threshs)

#%% 2. slope test: fold change in expected spots/pixel across the whole stack

def slope_fit(d):
    fit = smf.glm('n_spots ~ plane', data=d, family=sm.families.Poisson(),
                  offset=np.log(d['area'])).fit(
                      cov_type='cluster', cov_kwds={'groups': d['mask_id'].values})
    span = d['plane'].max() - d['plane'].min()
    b = fit.params['plane']
    lo, hi = fit.conf_int().loc['plane']
    return np.exp(b*span), np.exp(lo*span), np.exp(hi*span)   # fold change, 95% CI

print('all planes:      FC=%.3f  CI=(%.3f, %.3f)' % slope_fit(df))
print('without plane 0: FC=%.3f  CI=(%.3f, %.3f)' % slope_fit(df[df.plane >= 1]))
# FC is shallowest vs deepest. Bias is bounded if the CI sits inside, say, 0.9 to 1.1.

#%%
masks = pp_masks_im[row_slice, col_slice, :]
rows = []
for k in range(masks.shape[2]):
    m = masks[:,:,k]
    ids, areas = np.unique(m[m > 0], return_counts=True)
    rows += [(i, k, a) for i, a in zip(ids, areas)]
geo = pd.DataFrame(rows, columns=['Mask ID', 'plane', 'area'])
print('mask IDs in more than one plane:', geo['Mask ID'].duplicated().sum())  # expect 0
cells = cell_df.merge(geo, on='Mask ID')
span = cells.plane.max() - cells.plane.min()

def depth_fc(d, gene, extra=''):
    d = d.assign(y=d[gene])
    fit = smf.glm('y ~ plane' + extra, data=d, family=sm.families.Poisson(),
                  offset=np.log(d['area'])).fit(cov_type='HC1')
    b = fit.params['plane']; lo, hi = fit.conf_int().loc['plane']
    return np.exp(b*span), np.exp(lo*span), np.exp(hi*span)

for g in quant_genes:
    print('%-8s FC=%.3f  CI=(%.3f, %.3f)' % ((g,) + depth_fc(cells, g)))
    
#%%

import numpy as np, pandas as pd
import statsmodels.formula.api as smf
import statsmodels.api as sm
import matplotlib.pyplot as plt

g = 'meis2'
c = cells[cells['area'] > 0].copy()
c['y'] = c[g]
span = c.plane.max() - c.plane.min()
print(len(c), 'cells,', int(c.y.sum()), 'spots')

def fc(d, formula='y ~ plane', term='plane'):
    fit = smf.glm(formula, data=d, family=sm.families.Poisson(),
                  offset=np.log(d['area'])).fit(cov_type='HC1')
    b = fit.params[term]; lo, hi = fit.conf_int().loc[term]
    return np.exp(b*span), np.exp(lo*span), np.exp(hi*span)

#%% 1. linear slope, with and without plane 0
print('all planes:      FC=%.3f  CI=(%.3f, %.3f)' % fc(c))
print('without plane 0: FC=%.3f  CI=(%.3f, %.3f)' % fc(c[c.plane >= 1]))

#%% 2. shape: deep / middle / shallow thirds (catches a U that a slope misses)
c['third'] = pd.cut(c.plane, 3, labels=['deep', 'mid', 'shallow'])
fit = smf.glm('y ~ C(third, Treatment("mid"))', data=c, family=sm.families.Poisson(),
              offset=np.log(c['area'])).fit(cov_type='HC1')
ci = np.exp(fit.conf_int())
for name in fit.params.index[1:]:
    print('%-45s ratio vs mid = %.3f  CI=(%.3f, %.3f)' %
          (name, np.exp(fit.params[name]), ci.loc[name, 0], ci.loc[name, 1]))

#%% 3. per-plane density with bootstrap CI over cells
rng = np.random.default_rng(0)
planes = np.sort(c.plane.unique())
def dens(d):
    s = d.groupby('plane').agg(n=('y', 'sum'), a=('area', 'sum'))
    return (s.n / s.a).reindex(planes).values
boots = np.array([dens(c.iloc[rng.integers(0, len(c), len(c))]) for _ in range(1000)])
lo, hi = np.nanpercentile(boots, [2.5, 97.5], axis=0)
plt.figure(figsize=(7, 4))
plt.fill_between(planes, lo, hi, alpha=0.3)
plt.plot(planes, dens(c), 'o-')
plt.xlabel('plane (0 = deepest)'); plt.ylabel('in-mask spots per pixel')
plt.show()

#%% 4. power: inject a known trend by thinning spots, see how often it is detected
def inject(d, tc, deep_loses, rng):
    d = d.copy()
    frac = d.plane / d.plane.max()
    keep = 1 - tc * (1 - frac if deep_loses else frac)
    d['y'] = rng.binomial(d['y'].values, keep.values)
    return d

for tc in [0.10, 0.15, 0.20, 0.30]:
    for deep_loses in [True, False]:
        res = [fc(inject(c, tc, deep_loses, rng)) for _ in range(100)]
        hit = np.mean([(r[2] < 1) or (r[1] > 1) for r in res])
        print('%2.0f%% loss at %-7s end: detected in %3.0f%% of runs, mean FC=%.3f' %
              (tc*100, 'deep' if deep_loses else 'shallow', hit*100, np.mean([r[0] for r in res])))

#%%

def thirds(d, formula='y ~ C(third, Treatment("mid"))', offset=True):
    d = d[d.area > 0].copy()
    d['third'] = pd.cut(d.plane, 3, labels=['deep', 'mid', 'shallow'])
    fit = smf.glm(formula, data=d, family=sm.families.Poisson(),
                  offset=np.log(d['area']) if offset else None).fit(cov_type='HC1')
    ci = np.exp(fit.conf_int())
    for n in fit.params.index:
        if 'third' in n or 'area' in n:
            print('  %-42s %.3f  CI=(%.3f, %.3f)' % (n, np.exp(fit.params[n]), ci.loc[n,0], ci.loc[n,1]))
    return d

# 1. is spots ~ area really proportional? Free coefficient on log(area)
print('area free:')
d = thirds(c, 'y ~ C(third, Treatment("mid")) + np.log(area)', offset=False)
print(d.groupby('third', observed=True).area.median())

# 2. does the U depend on the threshold? (uses per_plane_table from earlier)
for s in [0.85, 1.0, 1.15]:
    t = per_plane_table(im, masks, fluo, threshs, scale=s).rename(columns={'n_spots': 'y'})
    print('threshold x%.2f:' % s); thirds(t)

# 3. does the U depend on the planes at the extremes? Drop 3 planes at each end
print('trimmed:'); thirds(c[(c.plane >= 3) & (c.plane <= c.plane.max() - 3)])

#%%

from bigfish.stack import log_filter
from bigfish.detection.spot_detection import local_maximum_detection

psf = img_utils.calc_psf_sigma(fluo) / 325
raw = quant_gene_ims[gene][row_slice, col_slice, :]
rows = []
for k in range(im.shape[2]):
    f = log_filter(im[:, :, k], psf)
    pk = local_maximum_detection(f, psf)
    med = np.median(f); mad = 1.4826 * np.median(np.abs(f - med))
    t = threshs[k]
    v_in = f[pk & (masks[:, :, k] > 0)]
    v_out = f[pk & (masks[:, :, k] == 0)]
    r = raw[:, :, k]
    rows.append(dict(plane=k,
                     thr_over_mad=(t - med) / mad,
                     n_in=(v_in > t).sum(), n_in_85=(v_in > .85*t).sum(),
                     n_in_115=(v_in > 1.15*t).sum(),
                     n_out=(v_out > t).sum(),
                     noise_gain=mad / (1.4826 * np.median(np.abs(r - np.median(r))) + 1e-9)))
diag = pd.DataFrame(rows)
diag['steepness'] = diag.n_in_85 / diag.n_in_115.clip(lower=1)
print(diag.round(2).to_string(index=False))

#%%

def u_ratios(d, trim=0, max_plane=37):
    d = d[(d.plane >= trim) & (d.plane <= max_plane - trim) & (d.area > 0)].copy()
    d['third'] = pd.cut(d.plane, 3, labels=['deep', 'mid', 'shallow'])
    fit = smf.glm('y ~ C(third, Treatment("mid")) + np.log(area)', data=d,
                  family=sm.families.Poisson()).fit(cov_type='HC1')
    ci = np.exp(fit.conf_int())
    out = []
    for n in fit.params.index[1:3]:
        out.append('%.2f (%.2f-%.2f)' % (np.exp(fit.params[n]), ci.loc[n, 0], ci.loc[n, 1]))
    return out

floors = []
for k in range(im.shape[2]):
    f = log_filter(im[:, :, k], psf)
    med = np.median(f); mad = 1.4826 * np.median(np.abs(f - med))
    floors.append((med, mad))
floors = np.array(floors)

def clamp(c):
    return np.maximum(threshs, floors[:, 0] + c * floors[:, 1])

# (planes 36-37 are at 7.05-8.6, so c=9 would also change the shallow end; stick to c <= 8)
for c in [7, 8]:
    t_c = clamp(c)
    for s in [0.85, 1.0, 1.15]:
        t = per_plane_table(im, masks, fluo, t_c, scale=s).rename(columns={'n_spots': 'y'})
        t = t[t.plane <= 37]
        dp, sh = u_ratios(t, 0)
        print('floor %d MAD  x%.2f  deep %s   shallow %s' % (c, s, dp, sh))
        
#%%

g = 'meis2'
cc = cells[cells['area'] > 0].copy()
cc['y'] = cc[g]
cc = cc[cc.plane <= 37]          # plane 38 has no in-mask data

pick = pd.concat([
    cc[cc.plane <= 8].sample(20, random_state=1),
    cc[(cc.plane >= 15) & (cc.plane <= 25)].sample(20, random_state=1)
])[['Mask ID', 'plane', 'area', 'y']]
print(pick.sort_values('plane').to_string(index=False))

#%%

from statsmodels.nonparametric.smoothers_lowess import lowess

whole = np.array(num_spots_dict['meis2'])[:38]          # drop plane 38
A_img = im.shape[0] * im.shape[1]
g = (cc.groupby('plane').agg(n_in=('y', 'sum'), a_in=('area', 'sum'))
       .reindex(range(38)).fillna(0))
g['whole'] = whole
# in-mask density relative to out-of-mask density; flat across plane => multiplicative bias
g['enrich'] = (g.n_in / g.a_in) / ((g.whole - g.n_in) / (A_img - g.a_in))
print(g.round(3).to_string())

# plane factor from the whole-image curve, relative to the middle of the stack
f = lowess(whole, np.arange(38), frac=0.3, return_sorted=False)
factor = f / np.median(f[9:30])
cc['y_adj'] = cc.y / factor[cc.plane.values]

# the check that matters: is the corrected in-mask count flat in plane?
d = cc.assign(y_raw=cc.y, y=cc.y_adj)
print('raw      ', u_ratios(cc.assign(y=cc.y_raw), 0))
print('corrected', u_ratios(d, 0))

#%%

import numpy as np, pandas as pd
import statsmodels.formula.api as smf
import statsmodels.api as sm
import matplotlib.pyplot as plt
from statsmodels.nonparametric.smoothers_lowess import lowess
from scipy.stats import spearmanr

N_PLANES = 38          # planes 0-37, plane 38 has no in-mask data

#%% helpers

def get_tables(gene, scale=1.0):
    """Per-cell table (mask_id, plane, area, y) and whole-image spot count per plane."""
    fl = gene_to_fluo[gene]
    img = quant_gene_ims[gene][row_slice, col_slice, :]
    img, _ = img_utils.correct_depth_attenuation(img, scale_to_log_filter=True,
                                                 plot_fit=False, fluo=fl)
    th = np.loadtxt(os.path.join(data_dir, gene + '_threshs.csv'), delimiter=',') * scale
    rows, whole = [], []
    for k in tqdm(range(N_PLANES)):
        spots = matrix_utils.count_spots(img[:, :, k], fl, thresh=th[k])
        whole.append(len(spots))
        m = masks[:, :, k]
        ids, areas = np.unique(m[m > 0], return_counts=True)
        if len(spots):
            hit = m[spots[:, 0], spots[:, 1]]
            cnt = pd.Series(hit[hit > 0]).value_counts()
        else:
            cnt = pd.Series(dtype=int)
        rows += [(i, k, a, int(cnt.get(i, 0))) for i, a in zip(ids, areas)]
    d = pd.DataFrame(rows, columns=['mask_id', 'plane', 'area', 'y'])
    return d, np.array(whole)

def correct(d, whole, frac=0.3, cap=False):
    """Divide each cell by a plane factor from the lowess-smoothed whole-image count."""
    f = lowess(whole, np.arange(len(whole)), frac=frac, return_sorted=False)
    factor = f / np.median(f[9:30])
    if cap:   # optional: never correct by more than the in-mask density implies
        s = d.groupby('plane').agg(n=('y', 'sum'), a=('area', 'sum')).reindex(range(N_PLANES))
        fm = lowess((s.n / s.a).fillna(0).values, np.arange(N_PLANES), frac=0.3, return_sorted=False)
        factor = np.minimum(factor, fm / np.median(fm[9:30]))
    return d.assign(y_adj=d.y / factor[d.plane.values]), factor

def slope_fc(d, col):
    d = d[d.area > 0].assign(yy=d[col])
    fit = smf.glm('yy ~ plane + np.log(area)', data=d,
                  family=sm.families.Poisson()).fit(cov_type='HC1')
    span = d.plane.max() - d.plane.min()
    lo, hi = fit.conf_int().loc['plane']
    return np.exp(fit.params['plane'] * span), np.exp(lo * span), np.exp(hi * span)

def ends_ratio(d, col):
    """Plane 0 and planes 36-37 vs the middle (planes 9-29), area exponent free."""
    grp = np.where(d.plane == 0, 'p0',
          np.where(d.plane >= 36, 'p36_37',
          np.where((d.plane >= 9) & (d.plane <= 29), 'mid', 'other')))
    d = d.assign(grp=grp, yy=d[col])
    d = d[(d.grp != 'other') & (d.area > 0)]
    fit = smf.glm('yy ~ C(grp, Treatment("mid")) + np.log(area)', data=d,
                  family=sm.families.Poisson()).fit(cov_type='HC1')
    ci = np.exp(fit.conf_int())
    return ['%s %.2f (%.2f-%.2f)' % (n.split('T.')[-1].rstrip(']'), np.exp(fit.params[n]),
                                     ci.loc[n, 0], ci.loc[n, 1])
            for n in fit.params.index if 'grp' in n]

def report(d, whole, label, frac=0.3, cap=False):
    d2, factor = correct(d, whole, frac, cap)
    print('--', label)
    for name, col in [('raw      ', 'y'), ('corrected', 'y_adj')]:
        dd = d2.assign(y=d2[col])
        print('  %s slope FC=%.3f (%.3f, %.3f) | thirds deep/shallow %s | ends %s' %
              ((name,) + slope_fc(d2, col) + (u_ratios(dd, 0), ends_ratio(d2, col))))
    return d2, factor

#%% 0. base case: meis2 at the original thresholds, plus cells per plane

d1, whole1 = get_tables('meis2', 1.0)
d1 = d1[d1.area > 0]

n_cells = d1.groupby('plane').size().reindex(range(N_PLANES), fill_value=0)
print(n_cells.to_string())
fit = smf.glm('n ~ plane', data=pd.DataFrame({'n': n_cells.values, 'plane': n_cells.index}),
              family=sm.families.Poisson()).fit(cov_type='HC1')
print('cells per plane: FC across stack=%.2f (%.2f, %.2f)' %
      (np.exp(fit.params['plane'] * 37), *np.exp(fit.conf_int().loc['plane'] * 37)))
print('spearman(plane, whole-image count):', spearmanr(np.arange(N_PLANES), whole1))

#%% 1 and 2. linear slope on corrected counts, and planes 0 / 36-37 directly

d1c, factor1 = report(d1, whole1, 'meis2 x1.00')
# if the ends are 10% or more below 1.00 above, rerun with the capped factor:
d1cap, _ = report(d1, whole1, 'meis2 x1.00, capped', cap=True)

# per-plane pooled density, normalised to the middle, with bootstrap CI over cells
rng = np.random.default_rng(0)
planes = np.arange(N_PLANES)
def dens(d, col):
    s = d.groupby('plane').agg(n=(col, 'sum'), a=('area', 'sum')).reindex(planes)
    return (s.n / s.a).values
def norm(v): return v / np.nanmedian(v[9:30])
boots = np.array([norm(dens(d1c.iloc[rng.integers(0, len(d1c), len(d1c))], 'y_adj'))
                  for _ in range(500)])
lo, hi = np.nanpercentile(boots, [2.5, 97.5], axis=0)
plt.figure(figsize=(8, 4))
plt.fill_between(planes, lo, hi, alpha=0.3, label='corrected, 95% CI')
plt.plot(planes, norm(dens(d1c, 'y_adj')), 'o-', label='corrected')
plt.plot(planes, norm(dens(d1c, 'y')), 'x--', label='raw')
plt.axhline(1, color='k', lw=0.5)
plt.xlabel('plane (0 = deepest)'); plt.ylabel('in-mask density / middle')
plt.legend(); plt.show()
# note: pooled spots per pixel is only a rough shape check, since spots scale with area^~2

#%% 3. sensitivity to the lowess smoothing
for fr in [0.2, 0.3, 0.5]:
    report(d1, whole1, 'meis2 x1.00, lowess frac=%.1f' % fr, frac=fr)

#%% 4. threshold robustness (about 1 minute per scale)
# the whole-image curve is recomputed at each scale so it moves with the cells
for s in [0.85, 1.15]:
    ds, ws = get_tables('meis2', s)
    ds = ds[ds.area > 0]
    report(ds, ws, 'meis2 x%.2f' % s)

#%% 5. second gene (needs <gene>_threshs.csv and the image loaded in quant_gene_ims)
GENE2 = 'syt1'     # change to whichever gene you quantify next
for s in [1.0, 0.85, 1.15]:
    d2g, w2g = get_tables(GENE2, s)
    d2g = d2g[d2g.area > 0]
    report(d2g, w2g, '%s x%.2f' % (GENE2, s))
    
#%%

all_masks = tifffile.imread(r'D:\2026-09-12_FN_SNr_2\post\masks\s01L\integrated_masks.tiff')
all_masks = np.transpose(all_masks,[1,2,0])
all_masks[np.isnan(all_masks)] = 0

#%% 1. imports and helpers
import os
import numpy as np, pandas as pd
import statsmodels.formula.api as smf
import statsmodels.api as sm
import matplotlib.pyplot as plt
from statsmodels.nonparametric.smoothers_lowess import lowess
from tqdm import tqdm

N_PLANES = 38                       # planes 0-37, plane 38 has no in-mask data
REF = slice(9, 30)                  # reference planes for normalisation (as before)
masks = pp_masks_im[row_slice, col_slice, :]               # curated cells
all_masks = all_masks[row_slice, col_slice, :]              # every segmented cell

def get_tables(gene, scale=1.0, exclude_curated=True):
    """One detection pass per plane. Returns the per-cell table, the whole-image count,
    and the spots-per-pixel density inside tissue (all cells; curated cells excluded
    by default so the factor is independent of the cells being corrected)."""
    fl = gene_to_fluo[gene]
    img = quant_gene_ims[gene][row_slice, col_slice, :]
    img, _ = img_utils.correct_depth_attenuation(img, scale_to_log_filter=True,
                                                 plot_fit=False, fluo=fl)
    th = np.loadtxt(os.path.join(data_dir, gene + '_threshs.csv'), delimiter=',') * scale
    rows, whole, tissue_n, tissue_a = [], [], [], []
    for k in tqdm(range(N_PLANES)):
        spots = matrix_utils.count_spots(img[:, :, k], fl, thresh=th[k])
        whole.append(len(spots))
        m = masks[:, :, k]
        ids, areas = np.unique(m[m > 0], return_counts=True)
        if len(spots):
            hit = m[spots[:, 0], spots[:, 1]]
            cnt = pd.Series(hit[hit > 0]).value_counts()
        else:
            cnt = pd.Series(dtype=int)
        rows += [(i, k, a, int(cnt.get(i, 0))) for i, a in zip(ids, areas)]

        tissue = all_masks[:, :, k] > 0
        if exclude_curated:
            tissue = tissue & (m == 0)
        tissue_a.append(tissue.sum())
        tissue_n.append(tissue[spots[:, 0], spots[:, 1]].sum() if len(spots) else 0)
    d = pd.DataFrame(rows, columns=['mask_id', 'plane', 'area', 'y'])
    d = d[d.area > 0]
    tissue_dens = np.array(tissue_n) / np.maximum(np.array(tissue_a), 1)
    return d, np.array(whole), tissue_dens

def norm_factor(v, frac=0.3):
    f = lowess(v, np.arange(len(v)), frac=frac, return_sorted=False)
    return f / np.median(f[REF])

def inmask_density(d):
    s = d.groupby('plane').agg(n=('y', 'sum'), a=('area', 'sum')).reindex(range(N_PLANES))
    return (s.n / s.a).fillna(0).values

def apply_factor(d, factor):
    return d.assign(y_adj=d.y / factor[d.plane.values])

def slope_fc(d, col):
    d = d[d.area > 0].assign(yy=d[col])
    fit = sm.GLM(d.yy, sm.add_constant(pd.DataFrame({'plane': d.plane, 'la': np.log(d.area)})),
                 family=sm.families.Poisson()).fit(cov_type='HC1')
    span = d.plane.max() - d.plane.min()
    lo, hi = fit.conf_int().loc['plane']
    return np.exp(fit.params['plane'] * span), np.exp(lo * span), np.exp(hi * span)

def u_ratios(d, col, trim=0, max_plane=37):
    d = d[(d.plane >= trim) & (d.plane <= max_plane - trim) & (d.area > 0)].copy()
    d['third'] = pd.cut(d.plane, 3, labels=['deep', 'mid', 'shallow'])
    d['yy'] = d[col]
    fit = smf.glm('yy ~ C(third, Treatment("mid")) + np.log(area)', data=d,
                  family=sm.families.Poisson()).fit(cov_type='HC1')
    ci = np.exp(fit.conf_int())
    return ['%.2f (%.2f-%.2f)' % (np.exp(fit.params[n]), ci.loc[n, 0], ci.loc[n, 1])
            for n in fit.params.index[1:3]]

def ends_ratio(d, col):
    """Deep end (planes 0-2), shallow end (planes 33-37) vs the middle (planes 9-29)."""
    grp = np.where(d.plane <= 2, 'deep_0_2',
          np.where(d.plane >= 33, 'shallow_33_37',
          np.where((d.plane >= 9) & (d.plane <= 29), 'mid', 'other')))
    d = d.assign(grp=grp, yy=d[col])
    d = d[(d.grp != 'other') & (d.area > 0)]
    fit = smf.glm('yy ~ C(grp, Treatment("mid")) + np.log(area)', data=d,
                  family=sm.families.Poisson()).fit(cov_type='HC1')
    ci = np.exp(fit.conf_int())
    return ['%s %.2f (%.2f-%.2f)' % (n.split('T.')[-1].rstrip(']'), np.exp(fit.params[n]),
                                     ci.loc[n, 0], ci.loc[n, 1])
            for n in fit.params.index if 'grp' in n]

#%% 2. does the all-cells factor track the in-mask factor? (about 1 minute)
GENE = 'meis2'
d1, whole1, tdens1 = get_tables(GENE, 1.0)

f_whole = norm_factor(whole1)
f_inmask = norm_factor(inmask_density(d1))
f_tissue = norm_factor(tdens1)

chk = pd.DataFrame({'whole': f_whole, 'tissue(all cells)': f_tissue, 'in_mask': f_inmask})
chk['implied_error_whole'] = chk['in_mask'] / chk['whole']
chk['implied_error_tissue'] = chk['in_mask'] / chk['tissue(all cells)']
print(chk.round(2).to_string())

fig, ax = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
ax[0].plot(f_whole, 'o-', label='whole image')
ax[0].plot(f_tissue, '^-', label='all cells (excl. curated)')
ax[0].plot(f_inmask, 's-', label='curated cells')
ax[0].axhline(1, color='k', lw=0.5); ax[0].set_ylabel('factor vs middle'); ax[0].legend()
ax[1].plot(chk.implied_error_whole, 'o-', label='whole image')
ax[1].plot(chk.implied_error_tissue, '^-', label='all cells')
ax[1].axhline(1, color='k', lw=0.5)
ax[1].set_ylabel('in-mask / factor'); ax[1].set_xlabel('plane (0 = deepest)'); ax[1].legend()
plt.tight_layout(); plt.show()

#%%

fw, fi = f_whole, f_inmask
def one_sided(fw, cap):
    f = fw.copy()
    k0 = int(np.argmin(fw))
    f[k0:] = np.minimum(f[k0:], cap)
    return f
for cap in [1.15, 1.2, 1.25, 1.3]:
    e = fi / one_sided(fw, cap)
    print('cap %.2f  planes 33-37: %s   worst: %.2f' % (cap, e[33:].round(2), np.abs(1 - e).max()))
#%%

from matplotlib.ticker import MaxNLocator

def plot_num_spots(num_spots, gene):
    z_plane = np.arange(len(num_spots)) + 1

    fig, ax = plt.subplots(figsize=(7, 4.5), dpi=120)

    ax.plot(
        z_plane, num_spots,
        color='#2E86AB', linewidth=2.2,
        marker='o', markersize=5,
        markerfacecolor='white', markeredgewidth=1.5,
        zorder=3
    )
    ax.fill_between(z_plane, num_spots, color='#2E86AB', alpha=0.08)

    ax.set_xlabel('Z-plane', fontsize=11, fontweight='medium')
    ax.set_ylabel('Number of spots detected', fontsize=11, fontweight='medium')
    ax.set_title(gene, fontsize=14, fontweight='bold', pad=12, style='italic')

    ax.xaxis.set_major_locator(MaxNLocator(integer=True, nbins=10))
    ax.grid(True, linestyle='--', alpha=0.4, zorder=0)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_ylim(bottom=0)

    fig.tight_layout()
    plt.show()
    
plot_num_spots(num_spots, 'tenm2')
    
    
    
#%%

cell_df.to_csv(exp_name + '_' + uniq_id + '_new_spot_all_cells' + '.csv',
               index=False)
    
#%%

import tifffile

def make_exp_img(pp_masks_im, cell_df, gene):
    
    max_mask = np.max(pp_masks_im)
    masks = cell_df['Mask ID'].values
    gene_exp = cell_df[gene].values
    
    lut = dict(zip(np.arange(max_mask), np.full(max_mask + 1, 0, dtype=np.uint16)))
    for i, mask in enumerate(masks):
        lut[mask] = np.uint16(gene_exp[i])
        
    array_lut = np.full(max_mask + 1, 0, dtype=np.uint16)
    array_lut[np.fromiter(lut.keys(), dtype=np.uint16)] = np.fromiter(lut.values(), dtype=np.uint16)
        
    exp_img = array_lut[pp_masks_im]
    exp_img = np.moveaxis(exp_img, 2, 0)
    
    tifffile.imwrite(
        gene + '.tif', 
        exp_img,
        metadata={'axes': 'ZYX'},
        imagej=True)
    
    
    
    
    
    
    
    
    



#%% threshold finding

i = 0
quant_gene = 'phox2b'
fluo = 'Alexa Fluor 594'

fluo = gene_to_fluo[quant_gene]
intensity_image = np.copy(quant_gene_ims[i])
intensity_image_crop = intensity_image[row_slice, col_slice, :]
intensity_image_crop = img_utils.correct_depth_attenuation(intensity_image_crop, 
                                                      scale_to_log_filter=True,
                                                      plot_fit=True,
                                                      fluo=fluo)

# thresholds_bigfish = matrix_utils.find_zstack_thresh_bigfish(intensity_image_crop, fluo)

#%%

from bigfish.detection.spot_detection import local_maximum_detection
import matplotlib.pyplot as plt

def mirror_peak_values(im_filtered, minimum_distance):
    vals = []
    for k in range(im_filtered.shape[2]):
        f = im_filtered[:,:,k]
        med = np.median(f)
        mirror = 2*med - f                      # reflect about the background median
        mask = local_maximum_detection(mirror, minimum_distance)
        vals.append(mirror[mask])
    return vals

def mirror_counts(vals, threshs):
    return np.array([(v > t).sum() if np.isfinite(t) else np.nan
                     for v, t in zip(vals, threshs)])

from scipy.stats import spearmanr

raw = np.array(first_thresh, dtype=float)
smooth = np.array(smooth_thresh, dtype=float)

mp = mirror_peak_values(im_filtered, minimum_distance)   # from the earlier code

tail99  = np.array([np.quantile(v, 0.99)  for v in mp])
tail999 = np.array([np.quantile(v, 0.999) for v in mp])

ok = np.isfinite(raw)
for name, s in [('q99', tail99), ('q999', tail999)]:
    print(name, 'vs raw knee   :', spearmanr(s[ok], raw[ok]))
    print(name, 'vs residual   :', spearmanr(s[ok], (raw - smooth)[ok]))
    
#%%

n_raw = mirror_counts(mp, raw)
n_smooth = mirror_counts(mp, smooth)
deep = np.arange(len(raw)) >= 2*len(raw)//3
print('deep planes, total mirror peaks  raw:', np.nansum(n_raw[deep]),
      ' smooth:', np.nansum(n_smooth[deep]))
print('shallow planes                   raw:', np.nansum(n_raw[~deep]),
      ' smooth:', np.nansum(n_smooth[~deep]))

# count at smooth vs raw threshold, using the real image, per plane
def real_counts(im_filtered, mask_fn, threshs):
    out = []
    for k, t in enumerate(threshs):
        f = im_filtered[:,:,k]
        out.append((f[mask_fn(f, minimum_distance)] > t).sum() if np.isfinite(t) else np.nan)
    return np.array(out)

c_raw = real_counts(im_filtered, local_maximum_detection, raw)
c_smooth = real_counts(im_filtered, local_maximum_detection, smooth)
print(np.c_[idx, c_raw, c_smooth, c_smooth / c_raw][::3])

#%%

import bigfish.detection as detection

from bigfish.plot import plot_detection

im = image[:,:,0]
fluo = 'Alexa Fluor 594'
thresh = 3517.9224

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
    

plot_detection(im, spots, rescale=True)

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
    
#%% threshold finding

i = 1
quant_gene = 'ralyl'
fluo = 'Alexa Fluor 647'

fluo = gene_to_fluo[quant_gene]
intensity_image = np.copy(quant_gene_ims[i])
intensity_image_crop = intensity_image[row_slice, col_slice, :]
intensity_image_crop = img_utils.correct_depth_attenuation(intensity_image_crop, 
                                                      scale_to_log_filter=True,
                                                      plot_fit=True,
                                                      fluo=fluo)
image = intensity_image_crop

#%%


    

    
        

    
    
#%%


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
