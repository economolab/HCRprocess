# -*- coding: utf-8 -*-
"""
Created on Fri Feb 14 15:18:09 2025

@author: jpv88
"""

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc

from tqdm import tqdm

from ABC_toolbox import cell_funcs

from sklearn.neighbors import LocalOutlierFactor

def find_nearest(array, value):
    
    array = array.astype(np.int64)
    idx = (np.abs(array - value)).argmin()
    
    return array[idx]


def bin_gene(gene, num_bins=10, thresh=2):
    """
    Bin gene expression values into equal-width (by magnitude) discrete intervals. 
    Removes outliers prior to binning. 

    Parameters
    ----------
    gene : numpy.ndarray (N)
        Gene expression vector of length N where N is the number of cells.
    num_bins : int, optional
        Total number of discrete bins to put expression values into.
        The default is 10.
    thresh : float, optional
        Number of standard deviations past which a point will be considered an outlier.
        The default is 3.
        
    Returns
    -------
    gene_out : numpy.ndarray (N)
        Gene expression vector of length N where N is the number of cells.
        Will now only be composed of (num_bins) total discrete values. 
    
    """
    
    if isinstance(gene, pd.Series):
        gene = gene.values
    
    clf = LocalOutlierFactor(metric='euclidean')
    clf.fit(gene.reshape(-1, 1))
    LOF = -1 * clf.negative_outlier_factor_
    
    # Determine which values exceed X times the standard deviation of the distribution.
    # zscore = stats.zscore(gene)
    # big_outlier = (zscore > thresh)
    # small_outlier = (zscore < -thresh)
    
    outlier = (LOF > thresh)

    # Bin only the non-outlier values
    out, bins = pd.cut(gene[~outlier], num_bins, retbins=True, labels=False)
    bins = (np.diff(bins) / 2) + bins[:-1]
    bins = np.round(bins)
    bins = bins.astype(np.uint16)
    out = bins[out]

    # Build output vector
    gene_out = gene.copy()
    gene_out[~outlier] = out
    
    outlier_vals = [find_nearest(bins, el) for el in gene[outlier]]
    gene_out[outlier] = np.array(outlier_vals)

    # Put the big outliers in the biggest bin and the small outliers in the smallest bin
    # gene_out[big_outlier] = bins[-1]
    # gene_out[small_outlier] = bins[0]

    return gene_out


def normalize_counts_to_median(exp):
    """
    Normalize all counts in a gene expression matrix to have the same total across cells.
    The target total is the median value across all cells.
    
    Parameters
    ----------
    exp : numpy.ndarray (n x m)
        Gene expression array with n cells and m genes.

    Returns
    -------
    exp_norm : numpy.ndarray (n x m)
        Gene expression array with n cells and m genes.
        All cells will now have the same total reads.

    """

    # calculate desired total reads per cell
    total_reads = np.sum(exp, axis=1, dtype=np.uint32)
    target_sum = np.median(total_reads)

    # normalize cells to target total reads
    adata = ad.AnnData(exp)
    sc.pp.normalize_total(adata, target_sum=target_sum, inplace=True)
    exp_norm = adata.X
    
    exp_norm = np.round(exp_norm)
    exp_norm = exp_norm.astype("uint16")

    return exp_norm


def calc_rdr(exp, meta, ratios, n=1000, N=10):
    """
    Calculate robust dynamic range of each gene.

    Parameters
    ----------
    exp : numpy.ndarray (n x m)
        Gene expression array with n cells and m genes.
    meta : pandas.DataFrame (n x m)
        Metadata DataFrame with n cells and m genes
    ratios : dict
        Relative proportions of each cluster in meta.
        Keys correspond to clusters.
        Values correspond to proportions (all sum to 1)
    n : int, optional
        Number of cells per bootstrapped distribution. 
        The default is 1000.
    N : int, optional
        Number of bootstrapped distributions to calculate RDR on. 
        The default is 10.

    Returns
    -------
    rdrs : numpy.ndarray (m)
        Vector of length m corresponding to RDR of each gene.

    """
    
    # list that stores RDRs for each gene from each bootstrapped distribution
    rdrs_total = []
    
    # iterate through bootstrapped distributions
    for i in tqdm(range(N)):

        boot_mat, boot_meta = cell_funcs.bootstrap_scRNAseq(
            meta, exp, ratios, n=n)
        
        # list that stores RDRs for each gene in this bootstrapped distribution
        rdrs_boot = []
        
        # iterate through genes
        for j in tqdm(range(boot_mat.shape[1])):
            
            gene_exp = boot_mat[:,j]

            gene_mean = np.mean(gene_exp)
            gene_std = np.std(gene_exp)
            
            # if outliers exist at top or bottom, cap max/min value to 3 stds off the mean
            # else set it to the biggest/smallest value 
            if np.any(gene_exp >= gene_mean + 3*gene_std):
                top = gene_mean + 3*gene_std
            else:
                top = np.max(gene_exp)

            if np.any(gene_exp <= gene_mean - 3*gene_std):
                bottom = gene_mean - 3*gene_std
            else:
                bottom = np.min(gene_exp)
            
            # calculate RDR for this gene
            rdrs_boot.append(top - bottom)

        rdrs_total.append(rdrs_boot)

    # take the mean across runs and log1p transform to smooth over outliers
    rdrs_total = np.array(rdrs_total)
    rdrs_total = np.mean(rdrs_total, axis=0)
    rdrs_total = np.log1p(rdrs_total)
    
    rdrs = rdrs_total

    return rdrs

def calc_var_E(exp, meta, ratios, n=1000, N=10, level='cluster'):
    
    var_E_total = []
    
    for i in tqdm(range(N)):
        
        boot_mat, boot_meta = cell_funcs.bootstrap_scRNAseq(
            meta, exp, ratios, n=n)

        meta_clu_dict = {}
        
        # build dictionary of expression arrays corresponding to each cluster 
        for clu in np.unique(meta[level]):
            meta_clu_dict[clu] = exp[(meta[level] == clu).values,:]
        
        frac_var_Es = []
        
        # iterate through genes
        for idx in tqdm(range(boot_mat.shape[1])):
            
            # dictionary that stores mean expression for each cluster for a given gene
            mean_clu_dict = {}
            
            # dictionary that stores variance of expression for each cluster for a given gene
            var_clu_dict = {}
            
            # iterate through clusters
            for clu in np.unique(meta[level]):
                
                mean_clu_dict[clu] = np.mean(meta_clu_dict[clu][:,idx])
                var_clu_dict[clu] = np.var(meta_clu_dict[clu][:,idx])
            
            var_UE = np.average(list(var_clu_dict.values()), weights=list(ratios.values()))
            var_E = np.var(list(mean_clu_dict.values()))
            
            frac_var_E = (var_E / (var_E + var_UE))
            
            if ~np.isnan(frac_var_E):
                frac_var_Es.append(frac_var_E)
            else:
                frac_var_Es.append(0)
                
        var_E_total.append(frac_var_Es)
        
    var_E_total = np.array(var_E_total)
    var_E_total = np.mean(var_E_total, axis=0)
    
    var_E = var_E_total
    
    return var_E

    
    
    
    
    
    
    