# -*- coding: utf-8 -*-
"""
Build convenience items that can be called upon in utils
"""

import os
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
import pandas as pd

# %%

def build_gene_MERFISH():
    
    f_gene_MERFISH = os.path.join(
        params.data_dir, "metadata\\MERFISH-C57BL6J-638850\\20231215"
    )
    
    f_gene_MERFISH = os.path.join(
        f_gene_MERFISH, "gene.csv"
    )
    
    gene = pd.read_csv(f_gene_MERFISH)
    
    gene = gene["gene_symbol"].values
    
    save_dir = os.path.join(dname, "util_files")
    file = os.path.join(save_dir, "gene_MERFISH.npy")
    np.save(file, gene)
    
def build_gene_scRNAseq():
    
    f_gene_scRNAseq = os.path.join(
        params.data_dir, "metadata\\WMB-10X\\20231215"
    )
    
    f_gene_scRNAseq = os.path.join(
        f_gene_scRNAseq, "gene.csv"
    )
    
    gene = pd.read_csv(f_gene_scRNAseq)
    
    gene = gene["gene_symbol"].values
    
    save_dir = os.path.join(dname, "util_files")
    file = os.path.join(save_dir, "gene_scRNAseq.npy")
    np.save(file, gene)
    
# %%
    
build_gene_MERFISH()
build_gene_scRNAseq()
    
    



