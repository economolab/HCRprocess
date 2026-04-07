# -*- coding: utf-8 -*-
"""
Created on Sat Mar 29 19:27:06 2025

@author: jpv88
"""

import pickle

import anndata as ad
import matplotlib as mpl
import matplotlib.pyplot as plt
import scanpy as sc

from ABC_toolbox import ABC_utils

def plot_gene(gene, meta, exp, level="cluster", filter_low=False, ratios=None):
    
    adata = ad.AnnData(exp)
    adata.var_names = ABC_utils.load_gene("scRNAseq")
    adata.var_names_make_unique()
    
    xlabels = list(meta[level].values)

    adata.obs["cluster"] = xlabels

    sc.pl.violin(adata, gene, groupby="cluster", rotation=-90, xlabel=level,
                 palette=ABC_utils.load_colors_dict(level))
    
    mpl.rcParams['image.composite_image'] = False
    plt.rcParams['svg.fonttype'] = 'none'
    
def fetch_colors_dict(level):

    with open("./ABC_toolbox/util_files/colors_dict.pkl", "rb") as handle:
        colors_dict = pickle.load(handle)

    match level:
        case "class":
            return_dict = colors_dict["class_dict"]
        case "subclass":
            return_dict = colors_dict["subclass_dict"]
        case "supertype":
            return_dict = colors_dict["supertype_dict"]
        case "cluster":
            return_dict = colors_dict["cluster_dict"]

    return return_dict
