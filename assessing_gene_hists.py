# -*- coding: utf-8 -*-
"""
Created on Sun Mar  2 19:06:56 2025

@author: jpv88
"""

from utils import utils, find_genes, classify_cells_2
from ABC_toolbox import cell_funcs, gene_funcs, iterative_reclustering
from sklearn.preprocessing import MinMaxScaler
import params
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import random
import pygad
import os
import sys

# obtain this file's directory, set it as the working directory, add it to the path
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)
os.chdir(dname)
sys.path.append(dname)

all_genes = utils.load_genes()

# %% load ant-IRN-PARN data

meta = pd.read_csv(os.path.join(params.util_files_dir,
                   "IRN_PARN_meta_pp.csv"), low_memory=False)
exp_raw = np.load(os.path.join(params.util_files_dir, "IRN_PARN_raw.npy"))
ratios = pd.read_pickle(os.path.join(
    params.util_files_dir, "IRN_PARN_ratios.pkl"))

# %% bootstrap super cells 

exp_super, meta_super, exp_super_raw = utils.boot_super(exp_raw, meta)

exp_super_raw = gene_funcs.normalize_counts_to_median(exp_super_raw)
exp_super_raw = np.round(exp_super_raw)
exp_super_raw = exp_super_raw.astype("uint16")

boot_mat, boot_meta = cell_funcs.bootstrap_scRNAseq(meta_super, 
                                                    exp_super_raw, 
                                                    ratios,
                                                    n=10000)

# %%

gene = 'Slc32a1'
gene_idx = np.where(all_genes == gene)[0][0]


# %% 

from sklearn.neighbors import LocalOutlierFactor

thresh = 2
clf = LocalOutlierFactor(metric='euclidean')
clf.fit(boot_mat[:,gene_idx].reshape(-1, 1))
LOF = -1 * clf.negative_outlier_factor_
outlier = (LOF > thresh)

scRNAseq_exp = boot_mat[:,gene_idx][~outlier]
scaler = MinMaxScaler()
scRNAseq_exp = scaler.fit_transform(scRNAseq_exp.reshape(-1, 1))

gene_df = pd.read_csv("s01-left_df.csv")

hcr_exp = gene_df[gene].values
clf = LocalOutlierFactor(metric='euclidean')
clf.fit(hcr_exp.reshape(-1, 1))
LOF = -1 * clf.negative_outlier_factor_
outlier = (LOF > thresh)

hcr_exp = hcr_exp[~outlier]
scaler = MinMaxScaler()
hcr_exp = scaler.fit_transform(hcr_exp.reshape(-1, 1))

# %%


plt.hist(scRNAseq_exp, bins=100)
plt.hist(hcr_exp, bins=100)
