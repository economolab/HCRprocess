# -*- coding: utf-8 -*-
"""
Created on Sat Feb 15 15:47:57 2025

@author: jpv88
"""

import numpy as np
import scanpy as sc

from sklearn.cluster import AgglomerativeClustering

from sklearn.metrics import accuracy_score
from sklearn.neighbors import KNeighborsClassifier
from sklearn.metrics import confusion_matrix
from sklearn.model_selection import StratifiedKFold

import random
import pandas as pd

from ABC_toolbox import classify_cells, cell_funcs

from anndata import AnnData
from tqdm import tqdm

import scipy as sp
from scipy.stats import entropy

# convert confusion matrix into distance matrix
# expects confusion matrix normalized over the true conditions (each row sums to 1)
def C2D(cm):
    
    cm_shape = cm.shape
    n_rows = cm_shape[0]
    n_cols = cm_shape[1]
    
    dm = np.zeros(cm_shape)
    
    # iterate through elements in the upper triangle of the confusion matrix
    for i in range(n_rows):
        for j in range(n_cols)[i+1:]:
            
            # each value is average confusion of a to b and b to a
            # inverted because greater confusion = less distance
            # mirrored across diagonal because distance is symmetric
            entry = 1 - ((cm[i,j] + cm[j,i]) / 2)
            dm[i,j] = entry
            dm[j,i] = entry
    
    return dm

# use a distance matrix 
def recluster(dm, step=0.001):
    
    n_labels = dm.shape[0]
    labels = range(n_labels)
    distance_threshold = step
    
    # keep going until you have less labels than you started with (at least 2 have been combined)
    while len(np.unique(labels)) == n_labels:
    
        clustering = AgglomerativeClustering(n_clusters=None,
                                             metric="precomputed",
                                             linkage="average",
                                             distance_threshold=distance_threshold)
    
        labels = clustering.fit_predict(dm)
        
        # on each iteration, increase the minimum distance at which two labels are combined
        distance_threshold += step    
    
    # determine which labels appear more than once (i.e. they represent combined labels)
    unique, unique_counts = np.unique(labels, return_counts=True)
    unique_multi = unique[unique_counts > 1]
    
    # build list of arrays where each array is indices of labels that have been combined
    new_groups = []
    for x in unique_multi:
        new_groups.append(np.array(range(n_labels))[labels == x])
    
    return new_groups

# manually build new_groups by inputting mix targets in a list of lists
def manual_recluster(labels, mix_targets):
    
    n_labels = len(labels)
    
    new_groups = []
    for mix in mix_targets:
        indices = np.isin(labels, mix)
        new_groups.append(np.array(range(n_labels))[indices])
        
    return new_groups

# builds list of lists where labels set to be combined by new_groups are in the same sublist 
def mix_labels(labels, new_groups):
    
    # convert labels list into a list of lists
    labels = list(labels)
    for i, val in enumerate(labels):
        temp = val
        labels[i] = []
        labels[i].append(temp)
    
    # labels that will be removed at the end, because they've been folded into another label
    idx_to_remove = []
    
    # iterate through new_groups appending their labels to the first label in the new_group
    # also mark the folded labels for deletion
    for group in new_groups:
        for i in range(len(group) - 1):
            labels[group[0]].append(labels[group[i+1]][0])
            idx_to_remove.append(group[i+1])
    
    # iterate backward through list of lists, deleting labels that have been folded into other labels
    for i in sorted(idx_to_remove, reverse=True):
        del labels[i]
    
    return labels

# returns a dict where keys are clusters and values are the cluster group that 
# cluster belongs to
def build_clu_mapping(labels, labels_mix):
    
    clu_mapping = {}
    
    # check if each label is in a mixed label group
    # if so, set it as a key in the clu_mapping with value equal to that mixed label group
    for label in labels:
        for label_list in labels_mix:
            if label in label_list:
                clu_mapping[label] = label_list
                break
    
    for k in clu_mapping.keys():
        clu_mapping[k] = str(clu_mapping[k])
        clu_mapping[k] = clu_mapping[k].translate({ord(c): None for c in "'[]"})
        
    keys = list(clu_mapping.keys())
    for k in keys:
        if ',' in k:
            k_split = k.split(", ")
            for k2 in k_split:
                clu_mapping[k2] = clu_mapping[k]
    
    return clu_mapping


def gini(array):
    """Calculate the Gini coefficient of a numpy array."""
    # based on bottom eq: http://www.statsdirect.com/help/content/image/stat0206_wmf.gif
    # from: http://www.statsdirect.com/help/default.htm#nonparametric_methods/gini.htm
    array = array.flatten() #all values are treated equally, arrays must be 1d
    if np.amin(array) < 0:
        array -= np.amin(array) #values cannot be negative
    array += 0.0000001 #values cannot be 0
    array = np.sort(array) #values must be sorted
    index = np.arange(1,array.shape[0]+1) #index per array element
    n = array.shape[0]#number of array elements
    return ((np.sum((2 * index - n  - 1) * array)) / (n * np.sum(array))) #Gini coefficient

    
# assumes raw counts from supercells
def gen_test_classifier(meta, exp, ratios, clu_mapping=None):
    
    if isinstance(exp, sp.sparse.sparray):
        exp = exp.todense()
        
    exp = pd.DataFrame(data=exp)
    
    meta = meta
    
    clusters = np.unique(meta["cluster"])
    
    def bootstrap_scRNAseq(meta, exp, ratios, n):
        
        meta.reset_index(drop=True, inplace=True)
        exp.reset_index(drop=True, inplace=True)
        
        num_per_clu = [round(x*n) for x in ratios.values()]

        boot_meta = []
        boot_mat = []
        
        # iterate through clusters, i is index, x is number of cells needed from that cluster
        for i, x in enumerate(num_per_clu):
            
            cur_clu = list(ratios.keys())[i]
            meta_clu = meta.iloc[meta.index[meta['cluster'] == cur_clu]]
            exp_clu = exp.iloc[exp.index[meta['cluster'] == cur_clu]]
            
            if len(meta_clu) != 0:
                for _ in range(x):
                    
                    idx = random.randint(0, len(meta_clu)-1)
                    boot_meta.append(meta_clu.iloc[idx])
                    boot_mat.append(exp_clu.iloc[idx])
                
        boot_meta = pd.DataFrame(boot_meta)
        boot_meta.reset_index(drop=True, inplace=True)
        boot_mat = pd.DataFrame(boot_mat)
        boot_mat.reset_index(drop=True, inplace=True)
        
        return boot_meta, boot_mat
    
    def prepare_bootstrapped_data(boot_meta, boot_mat):
        
        y = boot_meta["cluster"]
        
        X = boot_mat.to_numpy()
        
        return X, y
    
    def train_classifier(train, test, clu_mapping):
        
        exp_train = exp.iloc[train]
        exp_test = exp.iloc[test]
        
        meta_train = meta.iloc[train]
        meta_test = meta.iloc[test]
        
        meta_train, exp_train = bootstrap_scRNAseq(meta_train, exp_train, ratios, 4000)
        
        X_train, y_train = prepare_bootstrapped_data(meta_train, exp_train)
        
        if clu_mapping is not None:
            y_train = [clu_mapping[el] for el in y_train]
        
        # prepare test data
        meta_test, exp_test = bootstrap_scRNAseq(meta_test, exp_test, ratios, 1000)
        X_test, y_test = prepare_bootstrapped_data(meta_test, exp_test)
        
        if clu_mapping is not None:
            y_test = [clu_mapping[el] for el in y_test]
            
        combined_mat = np.vstack((X_train, X_test))
        adata = AnnData(combined_mat)
        adata.obs['batch'] = X_train.shape[0]*['train'] + X_test.shape[0]*['test']
        train_bool = (adata.obs['batch'] == 'train').values
        test_bool = (adata.obs['batch'] == 'test').values

        sc.pp.log1p(adata)
        sc.pp.scale(adata)
        sc.pp.pca(adata)
        
        pca_mat = adata.obsm['X_pca']
        train_pcs_mat = pca_mat[train_bool,:]
        test_pcs_mat = pca_mat[test_bool, :]
        
        neigh = KNeighborsClassifier(metric='cosine')
        neigh.fit(train_pcs_mat, y_train)
            
        final_pred = neigh.predict(test_pcs_mat)
        acc = accuracy_score(y_test, final_pred)
        
        labels = neigh.classes_
        cm = confusion_matrix(y_test, final_pred, labels=labels)
        cm = cm.astype(np.float64)
        
        new_ratios = {k: 0 for k in labels}
        for k in new_ratios.keys():
            k_split = k.split(", ")
            for k2 in k_split:
                new_ratios[k] += ratios[k2]
        
        ginis = []
        for i in range(cm.shape[0]):
            ginis.append(gini(cm[i,:]))
          
        sparsity = np.average(ginis, weights=list(new_ratios.values()))
        # sparsity = np.count_nonzero(cm==0) / np.size(cm)
        
        return acc, sparsity, cm, labels
    
    n_splits = 5
    skf = StratifiedKFold(n_splits=n_splits, shuffle=True)
    splits = skf.split(np.zeros(meta.shape[0]), meta['cluster'].values)
    indices = []
    for split in splits:
        indices.append(split[1])
    
    accs = []
    sparsities = []
    cms = []
    for i in tqdm(range(n_splits), desc='Train-test splits...'):
        train = indices[:i] + indices[i+1:]
        train = np.concatenate(train)
        train = np.sort(train)
        test = indices[i]
        
        acc, sparsity, cm, labels = train_classifier(train, test, clu_mapping)
        accs.append(acc)
        sparsities.append(sparsity)
        cms.append(cm)
    
    return np.mean(accs), np.mean(sparsities), cms, labels

def iterative_recluster_splits(boots, freqs, genes=None, acc_thresh=0.99, clu_mapping=None):
    
    final_cm = None
    final_labels = None
    
    all_accs = []
    all_n_labels = []
    bits = []
    
    acc = 0
    while acc < acc_thresh:
        
        full_cms = []
        accs = []
        
        data, cm_freqs = classify_cells.preprocess_data_splits(boots, freqs, clu_mapping=clu_mapping)
        
        res = classify_cells.cross_val_classifier_splits(data,  
                                                         cm_freqs,
                                                         genes=genes)
        labels = res['labels'][0]
        labels = [str(el) for el in labels]
        all_n_labels.append(len(labels))
        
        acc = np.mean(res["accs"])
        all_accs.append(acc)
        
        if clu_mapping is not None:
            cm_freqs = cell_funcs.freqs_to_cm_freqs(freqs, clu_mapping)
            S = entropy(list(cm_freqs.values()), base=2)
        else:
            S = entropy(list(freqs.values()), base=2)
        bits.append(S)
        
        cms = res['cms']
        full_cm = np.stack(cms)
        full_cm = np.sum(full_cm, axis=0)
        full_cms.append(full_cm)
        
        # convert confusion matrix to percentage per row
        full_cm = full_cm.astype(float)
        for i in range(full_cm.shape[0]):
            full_cm[i,:] = full_cm[i,:] / sum(full_cm[i,:])
        full_cm[np.isnan(full_cm)] = 0
        
        if acc >= acc_thresh:
            final_cm = full_cm
            final_labels = labels
            break
            
        dm = C2D(full_cm)
            
        new_groups = recluster(dm)
        new_labels = mix_labels(labels, new_groups)
        clu_mapping = build_clu_mapping(labels, new_labels)
        
        print(acc)
        print(len(new_labels))
        
        if len(new_labels) <= 1:
            final_cm = full_cm
            final_labels = labels
            break
    
    return final_labels, final_cm, all_accs, all_n_labels, bits

def iterative_recluster(meta, exp, ratios, acc_thresh=0.99):
    
    clu_mapping = None
    final_cm = None
    final_labels = None
    
    all_accs = []
    all_n_labels = []
    
    acc = 0
    while acc < acc_thresh:
        
        full_cms = []
        accs = []
        
        for i in range(5):
            acc, sparsity, cms, labels = gen_test_classifier(meta, 
                                                             exp, 
                                                             ratios,
                                                             clu_mapping=clu_mapping)
            
            labels = [str(el) for el in labels]
            accs.append(acc)
            full_cm = np.stack(cms)
            full_cm = np.sum(full_cm, axis=0)
            full_cms.append(full_cm)
        
        acc = np.mean(accs)
        all_accs.append(acc)
        all_n_labels.append(len(labels))
        
        full_cm = np.stack(full_cms)
        full_cm = np.sum(full_cm, axis=0)
        
        # convert confusion matrix to percentage per row
        full_cm = full_cm.astype(float)
        for i in range(full_cm.shape[0]):
            full_cm[i,:] = full_cm[i,:] / sum(full_cm[i,:])
        
        if acc >= acc_thresh:
            final_cm = full_cm
            final_labels = labels
            break
            
        dm = C2D(full_cm)
            
        new_groups = recluster(dm)
        new_labels = mix_labels(labels, new_groups)
        clu_mapping = build_clu_mapping(labels, new_labels)
        
        print(acc)
        print(len(new_labels))
        
        if len(new_labels) <= 1:
            final_cm = full_cm
            final_labels = labels
            break
    
    return final_labels, final_cm, all_accs, all_n_labels


