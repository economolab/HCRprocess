# -*- coding: utf-8 -*-
"""
Convenience functions
"""

import os
import pickle
import sys

# obtain this file's directory
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)

# add parent directory to path
pardir = os.path.dirname(dname)
sys.path.append(pardir)

# import params from parent directory
import params

import numpy as np

# %%

load_dir = os.path.join(dname, "util_files")

def load_gene(dtype):
    
    match dtype:
        
        case "MERFISH":
            file = os.path.join(load_dir, "gene_MERFISH.npy")
            gene = np.load(file, allow_pickle=True)
            return gene
        
        case "scRNAseq":
            file = file = os.path.join(load_dir, "gene_scRNAseq.npy")
            gene = np.load(file, allow_pickle=True)
            return gene
        
def load_parcel_dict():
    
    file = os.path.join(load_dir, "parcellation_dict.pickle")
    with open(file, 'rb') as handle:
        parc_dict = pickle.load(handle)

    return parc_dict
        
def load_colors_dict(level):

    with open(os.path.join(load_dir, "colors_dict.pkl"), "rb") as handle:
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

# identify corresponding taxonomic level of any taxonomic label
@np.vectorize
def id_tax(tax):

    with open("./ABC_toolbox/util_files/tax_dict.pkl", "rb") as handle:
        tax_dict = pickle.load(handle)

    for level in tax_dict.keys():
        if tax in tax_dict[level]:
            return level

    return "Unable to identify taxonomic label"

@np.vectorize
def fetch_cluster_tax(cluster):
    
    with open("./ABC_toolbox/util_files/cluster_dict.pkl", "rb") as handle:
        cluster_dict = pickle.load(handle)
    
    return cluster_dict[cluster] 


