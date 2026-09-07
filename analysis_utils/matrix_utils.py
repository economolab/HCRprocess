# -*- coding: utf-8 -*-
"""

"""

import bigfish.detection as detection
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import os

from analysis_utils import img_utils
from tqdm import tqdm
from sklearn.base import BaseEstimator, RegressorMixin
from sklearn.utils.validation import check_X_y, check_array, check_is_fitted
from sklearn.linear_model import RANSACRegressor
from sklearn.preprocessing import MinMaxScaler

import anndata as ad
import scanpy as sc

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


class QuadraticRegressor(BaseEstimator, RegressorMixin):
    """
    Fits a quadratic function y = a*x^2 + b*x + c using least squares.

    Works with a single feature (1D x). For multiple features it fits
    each feature's quadratic term plus a shared intercept via least squares
    on [x^2, x, 1].
    """

    def __init__(self):
        pass

    def fit(self, X, y):
        X, y = check_X_y(X, y, ensure_2d=True)
        if X.shape[1] != 1:
            raise ValueError(
                "QuadraticRegressor expects a single feature column (X.shape[1] == 1). "
                "Reshape your data with X.reshape(-1, 1) if needed."
            )

        x = X[:, 0]
        # Design matrix: [x^2, x, 1]
        A = np.column_stack([x**2, x, np.ones_like(x)])

        # Least squares solve
        coeffs, residuals, rank, singular_values = np.linalg.lstsq(A, y, rcond=None)

        self.a_, self.b_, self.c_ = coeffs
        self.coef_ = np.array([self.a_, self.b_])
        self.intercept_ = self.c_
        self.n_features_in_ = 1

        return self

    def predict(self, X):
        check_is_fitted(self, ["a_", "b_", "c_"])
        X = check_array(X, ensure_2d=True)
        x = X[:, 0]
        return self.a_ * x**2 + self.b_ * x + self.c_

    def score(self, X, y):
        # R^2 score (default from RegressorMixin uses this via sklearn's r2_score)
        from sklearn.metrics import r2_score
        return r2_score(y, self.predict(X))
    
def fit_threshold_curve(z_plane, threshs, plot_fit=False):
    
    reg = RANSACRegressor(estimator=QuadraticRegressor(),min_samples=0.5, max_trials=int(1e4))
    reg.fit(z_plane.reshape(-1, 1), threshs)
    smooth_thresh = reg.estimator_.predict(z_plane.reshape(-1,1))
    
    if plot_fit:
        
        plt.figure()
        plt.xlabel("Z-index")
        plt.ylabel("Spot counting threshold")
        plt.plot(z_plane, threshs, label='data')
        plt.scatter(z_plane[~reg.inlier_mask_], threshs[~reg.inlier_mask_], c='r', marker='x')
        plt.scatter(z_plane[reg.inlier_mask_], threshs[reg.inlier_mask_], c='g', marker='o')
        plt.plot(z_plane, reg.estimator_.predict(z_plane.reshape(-1, 1)), label='fitted')
        plt.legend()
    
    return smooth_thresh

def embed_cell_df(cell_df, n_components=2, embed="pca", 
                  ref_mask=None, norm_total=False, log_norm=False, z_score=True, 
                  clip_embedding=True, norm_embedding=True,
                  min_embed_pct=5, max_embed_pct=95):
    
    adata = ad.AnnData(cell_df.copy())
    
    if norm_total:
        sc.pp.normalize_total(adata)
    if log_norm:
        sc.pp.log1p(adata)
    if z_score:
        sc.pp.scale(adata)
    
    match embed:
        case "umap":
            obsm_key = "X_umap"
        case "pca":
            obsm_key = "X_pca"
        case _:
            raise ValueError(f"unknown embed: {embed!r}")
    
    if ref_mask is not None:
        
        adata_ref = adata[ref_mask].copy()
        adata_query = adata[~ref_mask].copy()
        
        match embed:
            case "umap":
                sc.pp.pca(adata_ref)
                sc.pp.neighbors(adata_ref)
                sc.tl.umap(adata_ref, n_components=n_components)
            case "pca":
                sc.pp.pca(adata_ref, n_comps=n_components)
                sc.pp.neighbors(adata_ref)
        
        sc.tl.ingest(adata_query, adata_ref, embedding_method=embed)
        
        adata_ref.obs["ref_or_query"] = "reference"
        adata_query.obs["ref_or_query"] = "query"
        adata = ad.concat(
            {"reference": adata_ref, "query": adata_query},
            label="ref_or_query",
            index_unique=None,
        )
        adata = adata[cell_df.index].copy()
        
    else:
            
        match embed:
            case "umap":
                sc.tl.pca(adata)
                sc.pp.neighbors(adata)
                sc.tl.umap(adata, n_components=n_components)
            case "pca":
                sc.tl.pca(adata, n_comps=n_components)
    
    cell_embedding = adata.obsm[obsm_key]
    
    if clip_embedding:
        embed_mins = np.percentile(cell_embedding, min_embed_pct, axis=0)
        embed_maxs = np.percentile(cell_embedding, max_embed_pct, axis=0)
        cell_embedding = np.clip(cell_embedding, a_min=embed_mins, a_max=embed_maxs)
        
    if norm_embedding:
        scaler = MinMaxScaler(feature_range=(0, 1))
        cell_embedding = scaler.fit_transform(cell_embedding)
    
    return cell_embedding


    
    
    
    


            
            