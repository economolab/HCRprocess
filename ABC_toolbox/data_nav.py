# -*- coding: utf-8 -*-
"""
Data navigation utility functions. 

@author: jpv88
"""

import h5py
import os
import pickle
import random
import scipy
import sys

import numpy as np
import pandas as pd
import nibabel as nib

from scipy.sparse import csr_array
from tqdm import tqdm

# obtain this file's directory
abspath = os.path.abspath(__file__)
dname = os.path.dirname(abspath)

# import params from the parent directory
pardir = os.path.dirname(dname)
sys.path.append(pardir)
import params


def load_meta(data_type):
    """
    Load in scRNAseq or MERFISH metadata

    Parameters
    ----------
    data_type: str, "scRNAseq" | "MERFISH"

    Returns
    -------
    meta: pd.DataFrame
        A pandas DataFrame of extracted scRNAseq or MERFISH metadata
        shape = (num_cells, num_columns)
    """

    if data_type == "scRNAseq":
        f_scRNAseq_meta = os.path.join(
            params.data_dir, "metadata\\WMB-10X\\20231215\\views"
        )
        f_scRNAseq_meta = os.path.join(
            f_scRNAseq_meta, "cell_metadata_with_cluster_annotation.csv"
        )

        print("Loading scRNAseq metadata...")

        meta = pd.read_csv(f_scRNAseq_meta, low_memory=False)

        print("Done!")

    elif data_type == "MERFISH":
        f_MERFISH = os.path.join(params.data_dir, "custom")
        f_MERFISH = os.path.join(f_MERFISH, "MERFISH_meta.csv")

        print("Loading MERFISH metadata...")

        meta = pd.read_csv(f_MERFISH, low_memory=False)

        print("Done!")

    return meta


def extract_MERFISH_meta(
    meta, selection, kind="restrict", category="anatomy", level="organ"
):
    """
    Extracts portions of a MERFISH metadata file

    Parameters
    ----------
    meta: pd.DataFrame
        A pandas DataFrame of MERFISH metadata
        shape = (num_cells, num_columns)
    selection: scalar or sequence
        Desired taxonomic types or anatomical regions
    kind: str, "restrict" | "remove", default "restrict"
        Whether to restrict the metadata to just the selection or remove the selection
    category: str, "anatomy" | "taxonomy", default "anatomy"
    level: str, if "anatomy", "organ" | "category" | "division" | "structure" | "substructure"
           str, if "taxonomy", "class" | "subclass" | "supertype" | "cluster"

    Returns
    -------
    meta: pd.DataFrame
        A pandas DataFrame of extracted MERFISH metadata
        shape = (num_cells, num_columns)
    """

    # convert to a list with one element if a scalar selection is input
    # (subsequent code expects a sequence)
    if np.isscalar(selection):
        selection = [selection]

    # dictionaries that convert input anatomy and taxonomy organization level
    # names into their actual metadata column names
    anatomy_level_dict = {
        "organ": "parcellation_organ",
        "category": "parcellation_category",
        "division": "parcellation_division",
        "structure": "parcellation_structure",
        "substructure": "parcellation_substructure",
    }

    tax_level_dict = {
        "class": "class",
        "subclass": "subclass",
        "supertype": "supertype",
        "cluster": "cluster",
    }

    # constructs a boolean mask corresponding to cells that match the input
    # selection in the given category at the given level
    def construct_bool_mask(meta, category, level):
        if category == "anatomy":
            bool_mask = np.zeros((len(meta)))

            for x in selection:
                bool_mask = bool_mask | (meta[anatomy_level_dict[level]] == x)

        elif category == "taxonomy":
            bool_mask = np.zeros((len(meta)))

            for x in selection:
                bool_mask = bool_mask | (meta[tax_level_dict[level]] == x)

        return bool_mask

    bool_mask = construct_bool_mask(meta, category, level)

    #
    if kind == "restrict":
        pass

    elif kind == "remove":
        bool_mask = np.logical_not(bool_mask)

    meta = meta.iloc[meta.index[bool_mask]]
    meta.reset_index(drop=True, inplace=True)

    return meta

def fetch_MERFISH(
    meta, form='raw'
):
    """
    Extracts portions of a MERFISH metadata file

    Parameters
    ----------
    meta: pd.DataFrame
        A pandas DataFrame of MERFISH metadata
        shape = (num_cells, num_columns)
    form: str, "log2" | "raw", default "raw"
        Expression matrix form, either log normalized or raw counts

    Returns
    -------
    exp: np.array
        A numpy array of expression data
        shape = (num_cells, num_genes)
    """
    
    MERFISH_data_dirs = [("MERFISH-C57BL6J-638850", "20230830", "C57BL6J-638850"),
                         ("Zhuang-ABCA-1", "20230830", "Zhuang-ABCA-1"),
                         ("Zhuang-ABCA-2", "20230830", "Zhuang-ABCA-2"),
                         ("Zhuang-ABCA-3", "20230830", "Zhuang-ABCA-3"),
                         ("Zhuang-ABCA-4", "20230830", "Zhuang-ABCA-4")]
    
    def build_f_MERFISH(data_dir, form):
        
        f_MERFISH_data = os.path.join(
            params.data_dir, "expression_matrices", data_dir[0], data_dir[1]
        )
        
        match form:
            case 'raw':
                f_MERFISH_data = os.path.join(
                    f_MERFISH_data, data_dir[2] + "-raw.h5ad"
                )
            case 'log2':
                f_MERFISH_data = os.path.join(
                    f_MERFISH_data, data_dir[2] + "-log2.h5ad"
                )
        
        return f_MERFISH_data
            
    # load expression matrix data and list of associated cell labels
    def load_exp_mat(f):
        
        f = h5py.File(f, "r")
        
        match type(f["X"]):
            
            case h5py._hl.group.Group:
        
                data = f["X"]["data"][:]
                indices = f["X"]["indices"][:]
                indptr = f["X"]["indptr"][:]
        
                exp_mat = csr_array((data, indices, indptr))
        
                cell_label = f["obs"]["cell_label"][:]
                cell_label = [str(x)[2:-1] for x in cell_label]
            
            # case h5py._hl.dataset.Dataset:
                
        return cell_label, exp_mat
    
    meta_cell_label = meta["cell_label"].values
    meta_list = []
    exp_list = []
    inds = []
    good_i = []
    i_to_find = list(range(len(meta_cell_label)))
    
    for data_dir in MERFISH_data_dirs:
        
        f_MERFISH_data = build_f_MERFISH(data_dir, form)
        
        print("Loading MERFISH data...")
        cell_label, exp_mat = load_exp_mat(f_MERFISH_data)
        
        cell_label = list(cell_label)
        cell_label = [str(x) for x in cell_label]
        
        for i in tqdm(i_to_find):
            try:
                ind = cell_label.index(meta_cell_label[i])
                inds.append(ind)
                meta_list.append(meta.iloc[i])
                good_i.append(i)
            except:
                pass
        
        exp = exp_mat[inds,:]
        exp_list.append(exp)
        
        for i in good_i:
            i_to_find.remove(i)
        
        inds = []
        good_i = []
        
    meta_full = pd.DataFrame(meta_list)
    
    return meta_full, exp


def fetch_scRNAseq(clusters, meta, form="log2", neur_frac=1, nn_frac=0.5):
    """
    Extracts all scRNAseq cells associated with any desired clusters

    Parameters
    ----------
    clusters: sequence
        Sequence of desired clusters
    meta: pd.DataFrame
        scRNAseq metadata
        shape = (num_cells, num_columns)
    form: str, "log2" | "raw", default "log2"
        Expression matrix form, either log normalized or raw counts

    Returns
    -------
    meta: pd.DataFrame
        A pandas DataFrame of extracted scRNAseq metadata
        shape = (num_cells, num_columns)
    exp: np.array
        A numpy array of expression data
        shape = (num_cells, num_genes)
    """

    # restrict metadata to input clusters
    bool_mask = np.zeros((len(meta)))
    
    print('Identifying which cells need to be found...')
    for cell_type in clusters:
        bool_mask = bool_mask | (meta["cluster"] == cell_type)
    print('Done!')
    meta = meta.iloc[meta.index[bool_mask]]
    meta.reset_index(drop=True, inplace=True)

    # calculates which expression matrix files need to be read based on input
    # sequence of possible cell rois
    def calc_files_to_read(rois):

        files_to_read = []
        
        for feat_mat in rois:
            files_to_read.append(feat_mat + '-' + form + '.h5ad')

        return files_to_read

    # load expression matrix data and list of associated cell labels
    def load_exp_mat(f):
        f = h5py.File(f, "r")

        data = f["X"]["data"][:]
        indices = f["X"]["indices"][:]
        indptr = f["X"]["indptr"][:]

        exp_mat = csr_array((data, indices, indptr))

        cell_label = f["obs"]["cell_label"][:]
        cell_label = [str(x)[2:-1] for x in cell_label]

        return cell_label, exp_mat
    
    def delete_row_csr(mat, i):

        n = mat.indptr[i+1] - mat.indptr[i]
        if n > 0:
            mat.data[mat.indptr[i]:-n] = mat.data[mat.indptr[i+1]:]
            mat.data = mat.data[:-n]
            mat.indices[mat.indptr[i]:-n] = mat.indices[mat.indptr[i+1]:]
            mat.indices = mat.indices[:-n]
        mat.indptr[i:-1] = mat.indptr[i+1:]
        mat.indptr[i:] -= n
        mat.indptr = mat.indptr[:-1]
        mat._shape = (mat._shape[0]-1, mat._shape[1])

    # take in current list of cell labels and expression matrix rows and add
    # any cells from a new file
    def load_file(meta_list, exp_mat_full, file, pbar, nn_frac):
        # set the appropriate path depending on scRNAseq data type
        if "10XMulti" in file:
            path = params.data_dir
            path = os.path.join(
                params.data_dir, "expression_matrices\\WMB-10XMulti\\20230830"
            )
        elif "10Xv2" in file:
            path = params.data_dir
            path = os.path.join(
                params.data_dir, "expression_matrices\\WMB-10Xv2\\20230630"
            )
        elif "10Xv3" in file:
            path = params.data_dir
            path = os.path.join(
                params.data_dir, "expression_matrices\\WMB-10Xv3\\20230630"
            )
            
        
        pbar.set_description('Loading matrix...' + ' (' + file + ')')
        cell_label, exp_mat = load_exp_mat(os.path.join(path, file))
        
        feat_mat_label = file
        feat_mat_label = feat_mat_label.replace('.h5ad', '')
        feat_mat_label = feat_mat_label.replace('-log2', '')
        feat_mat_label = feat_mat_label.replace('-raw', '')
        
        bool_mask = (meta['feature_matrix_label'] == feat_mat_label)
        
        # all the meta rows associated with cells in this file
        meta_file = meta.iloc[meta.index[bool_mask]]
        
        # non-neuronal classes
        non_neuronal = ['30 Astro-Epen', '31 OPC-Oligo', '32 OEC', '33 Vascular', '34 Immune']
        nn_mask = np.isin(meta_file["class"], non_neuronal)
        nn_index = list(meta_file.index[nn_mask])
        neur_index = list(meta_file.index[~nn_mask])
        
        # randomly keep some neurons
        k = round(neur_frac*len(neur_index))
        neur_keep = random.sample(neur_index, k)
        keep_index = list(neur_keep)
        
        # randomly keep some non-neuronal
        k = round(nn_frac*len(nn_index))
        nn_keep = random.sample(nn_index, k)
        keep_index.extend(nn_keep)
        
        meta_file = meta_file.loc[keep_index]
        
        # cell labels of cells of the desired clusters present in this matrix
        cells_present = meta_file['cell_label']
        
        pbar.set_description('Finding cells...(' + str(len(cells_present)) + ')')
        
        # indices in current file of cells that need to be extracted
        cells_present_indices = np.array([cell_label.index(x) for x in cells_present])
        exp_mat_full = scipy.sparse.vstack([exp_mat_full, exp_mat[cells_present_indices,:]])
        
        meta_list.append(meta_file)
        
        return meta_list, exp_mat_full

    # lists that contain cell labels and expression vectors as they're located
    meta_list = []

    # identify scRNAseq rois where input clusters are found
    rois = np.unique(meta["feature_matrix_label"])
    files_to_read = calc_files_to_read(rois)
    
    pbar = tqdm(
        total=len(files_to_read), position=0, leave=True
    )
    
    exp_mat_full = csr_array(np.ones((1, 32285), dtype='float32'))
    for file in files_to_read:
        cell_label_list, exp_mat_full = load_file(
                 meta_list, exp_mat_full, file, pbar, nn_frac
             )
        pbar.update()

    print("\nReformatting data...")
    
    exp_mat_full = exp_mat_full.tocsr()
    delete_row_csr(exp_mat_full, 0)

    meta_full = pd.concat(meta_list)

    print("Done!")

    return meta_full, exp_mat_full

def find_clusters(tax, level):

    file = open(os.path.join(dname, 'taxonomy_dict.pkl'), 'rb')
    taxonomy_dict = pickle.load(file)
    file.close()
    
    tax_level_dict = {
        "class": 0,
        "subclass": 1,
        "supertype": 2,
        "cluster": 3,
    }
    
    idx = tax_level_dict[level]
    
    keys = list(taxonomy_dict.keys())
    keys = np.array(keys)
    
    cluster_in_bool = np.zeros((1, len(keys)), dtype=bool)
    cluster_in_bool = np.squeeze(cluster_in_bool)

    clusters = []
    for i in range(len(cluster_in_bool)):
        if taxonomy_dict[keys[i]][idx] == tax:
            cluster_in_bool[i] = 1
            
    clusters = keys[cluster_in_bool]
    clusters = [str(x) for x in clusters]
    
    return clusters

def load_image_volume(vol="Allen"):
    
    if vol == "Allen":
        f_image_volume = os.path.join(
            params.data_dir, "image_volumes\\Allen-CCF-2020\\20230630"
        )
        f_image_volume = os.path.join(
            f_image_volume, "annotation_10.nii"
        )
        
        image_volume = nib.load(f_image_volume)
        image_volume = image_volume.get_fdata()
    
    elif vol == "MERFISH":
        f_image_volume = os.path.join(
            params.data_dir, "image_volumes\\MERFISH-C57BL6J-638850-CCF\\20230630"
        )
        f_image_volume = os.path.join(
            f_image_volume, "resampled_annotation.nii"
        )
        
        image_volume = nib.load(f_image_volume)
        image_volume = image_volume.get_fdata()
    
    return image_volume

def load_image_boundaries(vol="Allen"):
    
    if vol == "Allen":
        f_image_volume = os.path.join(
            params.data_dir, "image_volumes\\Allen-CCF-2020\\20230630"
        )
        f_image_volume = os.path.join(
            f_image_volume, "annotation_boundary_10.nii"
        )
        
        image_volume = nib.load(f_image_volume)
        image_volume = image_volume.get_fdata()
    
    elif vol == "MERFISH":
        f_image_volume = os.path.join(
            params.data_dir, "image_volumes\\MERFISH-C57BL6J-638850-CCF\\20230630"
        )
        f_image_volume = os.path.join(
            f_image_volume, "resampled_annotation_boundary.nii"
        )
        
        image_volume = nib.load(f_image_volume)
        image_volume = image_volume.get_fdata()
    
    return image_volume
    
    