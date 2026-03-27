# -*- coding: utf-8 -*-
"""
Take a fully completed HCRprocess image directory and turn it into a data
matrix.
"""

# %% imports

import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import tifffile

from tqdm import tqdm

#%% explicit and derived params

data_dir = 'D:\\2026-01-16_MC_SC_17\\post\\core_output\\s03L'
reg_tokens = ['nt', 'neur', 'snap25']
bin_tokens = ['bin']
mask_tokens = ['mask']
mask_im_tokens = ['masks_qc_final']

(head,uniq_id) = os.path.split(data_dir)
(post_dir,tail) = os.path.split(head)
(head,tail) =  os.path.split(post_dir)
(head,tail) = os.path.split(head)
exp_name = tail[11:]

#%% functions

# Check if any tokens are in a string. Returns bool
def check_for_tokens(string,tokens):
    
    TF = any([token.casefold() in string.casefold() for token in tokens])
    
    return TF

#%% load everything

data_f = np.array(os.listdir(data_dir))

reg_f_mask = [check_for_tokens(f,reg_tokens) for f in data_f]
reg_f = data_f[reg_f_mask]
data_f = data_f[np.invert(reg_f_mask)]

bin_f_mask = [check_for_tokens(f,bin_tokens) for f in data_f]
bin_f = data_f[bin_f_mask]
bin_f = str(bin_f[0])
data_f = data_f[np.invert(bin_f_mask)]

mask_f_mask = [check_for_tokens(f,mask_tokens) for f in data_f]
mask_f = data_f[mask_f_mask]
mask_f = str(mask_f[0])
data_f = data_f[np.invert(mask_f_mask)]

bin_df = pd.read_csv(os.path.join(data_dir,bin_f))
bin_df.drop(columns='Row',inplace=True)
bin_genes = list(bin_df.columns)

bin_genes_f_mask = [check_for_tokens(f,bin_genes) for f in data_f]
bin_genes_f = data_f[bin_genes_f_mask]
data_f = data_f[np.invert(bin_genes_f_mask)]

mask_df = pd.read_csv(os.path.join(data_dir,mask_f))

cell_df = pd.DataFrame()
cell_df.insert(0,'Cell ID',mask_df['Mask ID'])
cell_df.insert(0,'Mask ID',mask_df['Mask ID'])
cell_df.insert(cell_df.shape[1],'Z-plane',mask_df['Z-plane'])
cell_df.insert(cell_df.shape[1],'Z-span',mask_df['Z-span'])
cell_df.insert(cell_df.shape[1],'Assessed',mask_df['Assessed'])
cell_df.insert(cell_df.shape[1],'Result',mask_df['Result'])

for bin_gene in bin_genes:
    cell_df.insert(cell_df.shape[1],bin_gene,bin_df[bin_gene])
    cell_df[bin_gene] = cell_df[bin_gene].astype('bool')
    
cell_df = cell_df[cell_df['Assessed'] == 1]
cell_df = cell_df[cell_df['Result'] == 'good']
cell_df.drop(columns=['Assessed','Result'],inplace=True)
cell_df.reset_index(drop=True,inplace=True)

quant_genes = [f[3:-7] for f in data_f]

for quant_gene in quant_genes:
    cell_df.insert(cell_df.shape[1],quant_gene,np.zeros(cell_df.shape[0]))
    
cell_df['Cell ID'] = cell_df['Cell ID'].astype(str)
for i in range(cell_df.shape[0]):
    cell_df.loc[i,"Cell ID"] = exp_name + '_' + uniq_id + '_' + str(i+1)

# get the masks image path
mask_im_dir = os.path.join(post_dir,'masks',uniq_id)
mask_im_dir_f = np.array(os.listdir(mask_im_dir))
mask_im_f_mask = [check_for_tokens(f,mask_im_tokens) for f in mask_im_dir_f]
mask_im_f = mask_im_dir_f[mask_im_f_mask]
mask_im_f = str(mask_im_f[0])

    
#%% load images

quant_gene_ims = []

for i in tqdm(range(len(data_f)),desc='Loading images...'):
    im = tifffile.imread(os.path.join(data_dir,data_f[0]))
    im = np.transpose(im,[1,2,0])
    quant_gene_ims.append(im)
    
mask_im = tifffile.imread(os.path.join(mask_im_dir,mask_im_f))
mask_im = np.transpose(mask_im,[1,2,0])
    


