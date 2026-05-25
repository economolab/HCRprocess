# -*- coding: utf-8 -*-
"""
Utility functions for loading fully processed data for subsequent analysis
"""

import os

import anndata as ad
import numpy as np
import pandas as pd

def load_cell_df(filepath, usecols=None, index_col='Cell ID'):
    
    cell_df = pd.read_csv(filepath, index_col=index_col, usecols=usecols)
    
    for col in cell_df.columns:
        if not isinstance(cell_df[col].dtype, np.float64):
            cell_df[col] = cell_df[col].astype(np.float64)
    
    return cell_df

def concat_cell_dfs(cell_dfs):
    
    columns = [list(cell_df.columns) for cell_df in cell_dfs]
    if not all(cols == columns[0] for cols in columns):
        print("Columns do not match across all dataframes, concatenation failed.")
        return
    
    cell_df = pd.concat(cell_dfs)
    
    return cell_df

def load_experiment(exp_path, usecols=None, index_col='Cell ID'):
    
    (head, tail) = os.path.split(exp_path)
    exp_name = tail[11:]
    
    core_dir = os.path.join(exp_path, "post", "core_output")
    uniq_ids = os.listdir(core_dir)
    
    cell_dfs = []
    
    for uniq_id in uniq_ids:
        
        uniq_id_csv = exp_name + '_' + uniq_id + '.csv'
        filepath = os.path.join(core_dir, uniq_id, uniq_id_csv)
        
        if os.path.exists(filepath):
            cell_dfs.append(load_cell_df(filepath, usecols=usecols, index_col=index_col))
    
    return cell_dfs

def load_experiments(exp_paths, usecols=None, index_col='Cell ID'):
    
    cell_dfs = []
    
    for exp_path in exp_paths:
        
        cell_dfs.extend(load_experiment(exp_path, usecols=usecols, index_col=index_col))
        
    cell_df = concat_cell_dfs(cell_dfs)
    
    return cell_df
    
def build_cell_anndata(cell_df, obs_anno, norm_by_area=False):
    
    obs_anno = dict.fromkeys(obs_anno)
    
    for key in obs_anno.keys():
        obs_anno[key] = cell_df[key].values
        cell_df.drop(columns=key, inplace=True)
    
    if norm_by_area == True:
        
        one_pix = 0.325 * 0.325
        
        for col in cell_df.columns:
            cell_df[col] = cell_df[col] / (cell_df['Area'] * one_pix)
                                           
        cell_df.drop(columns='Area', inplace=True)
        
    adata = ad.AnnData(cell_df)
    
    adata.layers['raw'] = np.copy(adata.X)
    
    for key in obs_anno.keys():
        adata.obs[key] = obs_anno[key]
        
    return adata