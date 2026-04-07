# -*- coding: utf-8 -*-
"""
Utilities for classifying cell types
"""

import warnings

import numpy as np
import pandas as pd
import scanpy as sc
import scipy as sp

from anndata import AnnData, ImplicitModificationWarning
from multiprocessing import Pool, cpu_count
from sklearn.decomposition import PCA
from sklearn.metrics import accuracy_score, confusion_matrix
from sklearn.model_selection import StratifiedKFold
from sklearn.neighbors import KNeighborsClassifier
from sklearn.preprocessing import StandardScaler
from tqdm import tqdm

from ABC_toolbox import ABC_utils, cell_funcs

# %%

# train a cell type classifier
def train_classifier(X_train, y_train, clf_method='knn'):
    

    match clf_method:
        case 'knn':
            clf = KNeighborsClassifier(metric='cosine')

                
    clf.fit(X_train, y_train)

    return clf

def train_classifier_single(X_train, y_train, clf_method='knn'):
    
    match clf_method:
        case 'knn':
            clf = KNeighborsClassifier(metric='euclidean')

    clf.fit(X_train, y_train)

    return clf

# test a cell type classifier
def test_classifier(X_test, y_test, clf, freqs):
        
    y_pred = clf.predict(X_test)
    
    acc = accuracy_score(y_test, y_pred)
    
    labels = clf.classes_
    cm = confusion_matrix(y_test, y_pred, labels=labels)
    cm = cm.astype(np.float64)
    
    ginis = []
    for i in range(cm.shape[0]):
        ginis.append(gini(cm[i, :]))
    
    sparsity = np.average(ginis, weights=[freqs[label] for label in labels])
    
    return acc, sparsity, cm, labels
    

# make sure bootstrapped data is ready to be input to a classifier
def preprocess_data(exp, meta, freqs, clu_mapping=None):
    
    if isinstance(exp, sp.sparse.sparray):
        exp = exp.todense()
        
    # convert freqs to cluster mapping specific freqs
    if clu_mapping is not None:
        freqs = cell_funcs.freqs_to_cm_freqs(freqs, clu_mapping)

    y = meta["cluster"]

    if not isinstance(exp, np.ndarray):
        exp = exp.to_numpy()
    else:
        pass 
        
    exp = pd.DataFrame(data=exp)

    adata = AnnData(exp)

    sc.pp.log1p(adata)
    sc.pp.scale(adata)
    
    all_genes = ABC_utils.load_gene("scRNAseq")
    adata.var_names = all_genes
    adata.var_names_make_unique()
    
    X = adata

    return X, y, freqs

# preprocess expression matrix
def preprocess_exp(exp, scaler=None, return_scaler=False):
    
    # check if scipy sparse array
    if isinstance(exp, sp.sparse.sparray):
        exp = exp.todense()
    
    # make sure exp is a numpy array
    if not isinstance(exp, np.ndarray):
        exp = exp.to_numpy()
    else:
        pass 
    
    # log1p transform
    exp = np.log1p(exp)
    
    # make new scaler if no scaler is provided, otherwise use provided scaler
    if scaler == None:
        scaler = StandardScaler()
        exp = scaler.fit_transform(exp)
    else: 
        exp = scaler.transform(exp)
    
    # return scaler exp and scaler if requested, otherwise just exp
    if return_scaler == True:
        return exp, scaler
    else:
        return exp
    

# make sure bootstrapped data is ready to be input to a classifier
def preprocess_data_splits(boots, freqs, clu_mapping=None):
    
    # convert freqs to cluster mapping specific freqs
    if clu_mapping is not None:
        freqs = cell_funcs.freqs_to_cm_freqs(freqs, clu_mapping)
    
    data = []
    
    for boot in tqdm(boots, desc="Preprocessing bootstrapped data..."):
    
        train_exp = boot[0][0]
        train_meta = boot[0][1]
        test_exp = boot[1][0]
        test_meta = boot[1][1]
        
        train_exp, train_scaler = preprocess_exp(train_exp, return_scaler=True)
        test_exp = preprocess_exp(test_exp, scaler=train_scaler)
        
        exp = np.concatenate((train_exp, test_exp), axis=0)

        batch = np.concatenate((np.repeat('train', len(train_meta)),
                                np.repeat('test', len(test_meta))))
        
        meta = pd.concat((train_meta, test_meta))
        meta['batch'] = batch
        
        if clu_mapping is not None:
            new_clus = [clu_mapping[clu] for clu in meta['cluster'].values]
            meta['cluster'] = new_clus
        
        data.append((exp, meta))
        
    return data, freqs


def gini(array):
    """Calculate the Gini coefficient of a numpy array."""

    # from: http://www.statsdirect.com/help/default.htm#nonparametric_methods/gini.htm

    array = array.flatten()  # all values are treated equally, arrays must be 1d

    if np.amin(array) < 0:
        array -= np.amin(array)  # values cannot be negative

    array += 0.0000001  # values cannot be 0
    array = np.sort(array)  # values must be sorted

    index = np.arange(1, array.shape[0]+1)  # index per array element
    n = array.shape[0]  # number of array elements

    # Gini coefficient
    return ((np.sum((2 * index - n - 1) * array)) / (n * np.sum(array)))

# cross-validate a classifier
def cross_val_classifier_splits(data, freqs,
                                genes=None, clf_method='knn', verbose=True):
    
    """
    Parameters
    ----------
    meta : TYPE
        DESCRIPTION.
    exp : TYPE
        DESCRIPTION.
    freqs : TYPE
        DESCRIPTION.
    genes : TYPE, optional
        DESCRIPTION. The default is None.
    clf_method : TYPE, optional
        DESCRIPTION. The default is 'knn'.
    clu_mapping : TYPE, optional
        DESCRIPTION. The default is None.
    n_splits : TYPE, optional
        DESCRIPTION. The default is 5.

    Returns
    -------
    None.

    """
    
    n_splits = len(data)
    
    split_dicts = []
        
    for i in tqdm(range(n_splits), desc="Preprocessing splits...") if verbose else range(n_splits):
        
        exp = data[i][0]
        meta = data[i][1]
        
        train_mask = (meta['batch'].values == 'train')
        test_mask = (meta['batch'].values == 'test')
        gene_mask = np.array([gene in genes for gene in ABC_utils.load_gene("scRNAseq")])
        
        exp_train = exp[train_mask,:]
        exp_test = exp[test_mask,:]
        
        exp_train = exp_train[:,gene_mask]
        exp_test = exp_test[:,gene_mask]
        
        n_components = min(exp_train.shape) - 1
        if n_components > 50:
            n_components = 50
        
        pca = PCA(n_components=n_components)
        exp_train = pca.fit_transform(exp_train)
        exp_test = pca.transform(exp_test)
        
        split_dict = {}
        split_dict['X_train'] = exp_train
        split_dict['X_test'] = exp_test
        split_dict['y_train'] = meta['cluster'].values[train_mask]
        split_dict['y_test'] = meta['cluster'].values[test_mask]
        
        split_dicts.append(split_dict)
    
    res = {}
    res["accs"] = []
    res["sparsities"] = []
    res["cms"] = []
    res["labels"] = []
    
    for i in tqdm(range(n_splits), desc='Train-test splits...') if verbose else range(n_splits):
        split_dict = split_dicts[i]
        
        X_train = split_dict["X_train"]
        y_train = split_dict["y_train"]
        X_test = split_dict["X_test"]
        y_test = split_dict["y_test"]
        
        clf = train_classifier(X_train, y_train, 
                               clf_method=clf_method)
        
        acc, sparsity, cm, labels = test_classifier(X_test, y_test, clf, freqs)
        
        res["accs"].append(acc)
        res["sparsities"].append(sparsity)
        res["cms"].append(cm)
        res["labels"].append(labels)
        
    return res

# cross-validate a classifier
def cross_val_classifier(meta, exp, freqs,
                         genes=None, clf_method='knn', n_splits=5,
                         n_cells_boot=5000, verbose=True):
    
    """
    Parameters
    ----------
    meta : TYPE
        DESCRIPTION.
    exp : TYPE
        DESCRIPTION.
    freqs : TYPE
        DESCRIPTION.
    genes : TYPE, optional
        DESCRIPTION. The default is None.
    clf_method : TYPE, optional
        DESCRIPTION. The default is 'knn'.
    clu_mapping : TYPE, optional
        DESCRIPTION. The default is None.
    n_splits : TYPE, optional
        DESCRIPTION. The default is 5.

    Returns
    -------
    None.

    """
    
    # restrict to only a subset of genes
    if genes is not None:
        genes = np.array(genes)
        exp = exp[:,genes]
    
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message="Setting element")
        
        if len(genes) > 2:
            sc.pp.pca(exp)
    
    if len(genes) > 2:
        pcs = exp.obsm['X_pca']
    else:
        pcs = exp.X
   
    skf = StratifiedKFold(n_splits=n_splits, shuffle=True)
    splits = skf.split(np.zeros(meta.shape[0]), meta)
    indices = []
    for split in splits:
        indices.append(split[1])
    
    res = {}
    res["accs"] = []
    res["sparsities"] = []
    res["cms"] = []
    res["labels"] = []
    
    for i in (tqdm(range(n_splits), desc='Train-test splits...') if verbose else range(n_splits)):
        train = indices[:i] + indices[i+1:]
        train = np.concatenate(train)
        train = np.sort(train)
        test = indices[i]
        
        X_train = pcs[train,:]
        y_train = meta.iloc[train]
        X_test = pcs[test,:]
        y_test = meta.iloc[test]
        
        clf = train_classifier(X_train, y_train, 
                               clf_method=clf_method)
        
        acc, sparsity, cm, labels = test_classifier(X_test, y_test, clf, freqs)
        
        res["accs"].append(acc)
        res["sparsities"].append(sparsity)
        res["cms"].append(cm)
        res["labels"].append(labels)
        
    return res

# cross-validate a classifier
def cross_val_classifier_single(meta, exp, freqs,
                         gene, clf_method='knn', n_splits=5,
                         n_cells_boot=5000, verbose=True):
    
    """
    Parameters
    ----------
    meta : TYPE
        DESCRIPTION.
    exp : TYPE
        DESCRIPTION.
    freqs : TYPE
        DESCRIPTION.
    genes : TYPE, optional
        DESCRIPTION. The default is None.
    clf_method : TYPE, optional
        DESCRIPTION. The default is 'knn'.
    clu_mapping : TYPE, optional
        DESCRIPTION. The default is None.
    n_splits : TYPE, optional
        DESCRIPTION. The default is 5.

    Returns
    -------
    None.

    """
    
    exp = exp[:,gene]
    
    pcs = exp.X
   
    skf = StratifiedKFold(n_splits=n_splits, shuffle=True)
    splits = skf.split(np.zeros(meta.shape[0]), meta)
    indices = []
    for split in splits:
        indices.append(split[1])
    
    res = {}
    res["accs"] = []
    res["sparsities"] = []
    res["cms"] = []
    res["labels"] = []
    
    for i in (tqdm(range(n_splits), desc='Train-test splits...') if verbose else range(n_splits)):
        train = indices[:i] + indices[i+1:]
        train = np.concatenate(train)
        train = np.sort(train)
        test = indices[i]
        
        X_train = pcs[train,:]
        y_train = meta.iloc[train]
        X_test = pcs[test,:]
        y_test = meta.iloc[test]
        
        clf = train_classifier_single(X_train, y_train, 
                               clf_method=clf_method)
        
        acc, sparsity, cm, labels = test_classifier(X_test, y_test, clf, freqs)
        
        res["accs"].append(acc)
        res["sparsities"].append(sparsity)
        res["cms"].append(cm)
        res["labels"].append(labels)
        
    return res
