# -*- coding: utf-8 -*-
"""
Created on Sun Feb 23 15:05:22 2025

@author: jpv88
"""

import random
import numpy as np
import pandas as pd

from ABC_toolbox import gene_funcs

from sklearn.model_selection import StratifiedKFold
from tqdm import tqdm

def bootstrap_scRNAseq(meta, exp, freqs, n=1000):
    
    """
    Bootstrap a new distribution of scRNA-seq data where proportions match prescribed values.
    These proportions are calculated using MERFISH data.
    
    Parameters
    ----------
    meta : pandas.DataFrame (l x m)
        Metadata DataFrame with l cells and m genes
    exp : numpy.ndarray (l x m)
        Expression array with l cells and m genes
    ratios : dict
        Relative proportions of each cluster in meta.
        Keys correspond to clusters.
        Values correspond to proportions (all sum to 1)
    n : int
        Total number of cells desired in bootstrapped distribution.

    Returns
    -------
    boot_mat : numpy.ndarray (n x m)
        Bootstrapped expression array with n cells and m genes.

    """
    
    # buggy indices can mess everything up, reset to avoid
    meta.reset_index(drop=True, inplace=True)

    # determine how many cells of each cluster are needed
    num_per_clu = [round(x*n) for x in freqs.values()]
    
    # lists that stores 1D numpy arrays corresponding to each bootstrapped cell
    boot_mat = []
    boot_meta = []

    # iterate through clusters, i is index, x is number of cells needed from that cluster
    for i, x in enumerate(num_per_clu):

        # extract metadata and expression values associated with a given cluster
        cur_clu = list(freqs.keys())[i]
        meta_clu = meta.iloc[meta.index[meta['cluster'] == cur_clu]]
        exp_clu = exp[meta['cluster'] == cur_clu,:]
        
        # make sure the cluster exists
        if len(meta_clu) != 0:
            for _ in range(x):
                idx = random.randint(0, len(meta_clu)-1)
                boot_mat.append(exp_clu[idx,:])
                boot_meta.append(meta_clu.iloc[idx])

    boot_mat = np.array(boot_mat)
    boot_meta = pd.DataFrame(boot_meta)
    boot_meta.reset_index(drop=True, inplace=True)

    return boot_mat, boot_meta

def bootstrap_scRNAseq_splits(supers, freqs, n=1000):
    
    """

    """
    
    n_splits = len(supers)
    n_train = round(((n_splits - 1) / n_splits) * n)
    n_test = round((1 / n_splits) * n)
    
    boots = []
    
    freqs_cutoff = n_splits / n
    freqs = {k: v for k, v in freqs.items() if v >= freqs_cutoff}
    freqs = {k: v/sum(freqs.values()) for k, v in freqs.items()}
    
    for super_ in tqdm(supers, desc="Bootstrapping data..."):
        
        exp_train = super_[0][0]
        meta_train = super_[0][1]
        exp_test = super_[1][0]
        meta_test = super_[1][1]
        
        exp_train_boot, meta_train_boot = bootstrap_scRNAseq(meta_train, 
                                                             exp_train, 
                                                             freqs, 
                                                             n=n_train)
        
        exp_test_boot, meta_test_boot = bootstrap_scRNAseq(meta_test, 
                                                           exp_test, 
                                                           freqs, 
                                                           n=n_test)
        
        train_tuple = (exp_train_boot, meta_train_boot)
        test_tuple = (exp_test_boot, meta_test_boot)
        
        boots.append((train_tuple, test_tuple))
        
    return boots

def calc_frac_per_type(meta, level='cluster'):
    
    # find every unique cell type of the prescribed level in this metadata dataframe
    uniq_clu = np.unique(meta[level])
    
    # find how many cells of each type there are
    num_per_clu = []
    for clu in uniq_clu:
        num_per_clu.append(sum(meta[level] == clu))
    
    # calculate the fraction per cell type
    num_cells = len(meta)
    frac_per_clu = [x/num_cells for x in num_per_clu]
    
    # convert to dictionary
    d = {uniq_clu[i]: frac_per_clu[i] for i in range(len(uniq_clu))}
    
    return d

# convert MERFISH frequencies to new frequencies based on a cluster mapping
def freqs_to_cm_freqs(freqs, clu_mapping):
    
    cgs = np.unique(list(clu_mapping.values()))
    
    cm_freqs = {cg: 0 for cg in cgs}

    for cg in cm_freqs.keys():
        cg_split = cg.split(", ")
        for clu in cg_split:
            cm_freqs[cg] += freqs[clu]
            
    return cm_freqs

# bootstrap supercells
def boot_super(exp, meta, k=5, boot_factor=2):
    
    exp = exp.astype("uint16")
    
    cluster = meta["cluster"].values
    uniq_clu = np.unique(cluster)
    
    exp_super = []
    meta_super = []
    
    for clu in tqdm(uniq_clu, desc="Bootstrapping supercells for each cluster..."):
        
        mask = (cluster == clu)
        meta_clu = meta[mask]
        meta_clu.reset_index(inplace=True)
        exp_clu = exp[mask,:]
        n_real_cells = exp_clu.shape[0]
        
        for _ in range(max(1, round(n_real_cells*boot_factor))):
            
            idx = random.choices(range(n_real_cells), k=k)
            exps = [exp_clu[i,:] for i in idx]
            exp_super.append(sum(exps))
            meta_super.append(meta_clu.iloc[idx[0]])
            
    exp_super = np.vstack(exp_super)
    meta_super = pd.DataFrame(meta_super)
    
    return exp_super, meta_super

# boot super cells with pre split data
def boot_super_splits(exp, meta, splits, k=5, boot_factor=2):
    
    supers = []
    
    for split in tqdm(splits, desc="Boostrapping supercells for a split..."):
        
        train = split[0]
        test = split[1]
        
        exp_train = exp[train,:]
        meta_train = meta.iloc[train]
        exp_test = exp[test,:]
        meta_test = meta.iloc[test,:]
        
        exp_super_train, meta_super_train = boot_super(exp_train, 
                                                       meta_train, 
                                                       k=k, 
                                                       boot_factor=boot_factor)
        
        exp_super_test, meta_super_test = boot_super(exp_test, 
                                                     meta_test, 
                                                     k=k, 
                                                     boot_factor=boot_factor)
        
        train_tuple = (exp_super_train, meta_super_train)
        test_tuple = (exp_super_test, meta_super_test)
        
        supers.append((train_tuple, test_tuple))
    
    return supers

# filter cell types that don't have at least k cells
def filt_cells(exp, meta, k=5):
    
    meta.reset_index(drop=True, inplace=True)
    
    clu, clu_count = np.unique(meta["cluster"], return_counts=True)
    
    bool_mask = (clu_count >= k)
    keep_clu = clu[bool_mask]
    bool_mask = np.isin(meta["cluster"], keep_clu)
    
    meta = meta.iloc[meta.index[bool_mask]]
    meta.reset_index(drop=True, inplace=True)
    
    exp = exp[bool_mask,:]
    
    return exp, meta

# K fold train-test split cells
def K_fold_cells(meta, k=5, level='cluster'):
    
    skf = StratifiedKFold(n_splits=k, shuffle=True)
    
    # total number of cells
    n_samples = len(meta)
    
    X = np.zeros(n_samples)
    y = meta[level].values
    
    splits_gen = skf.split(X, y)
    
    splits = []
    for split in splits_gen:
        splits.append(split)
    
    return splits    
    
    
