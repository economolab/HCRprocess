# -*- coding: utf-8 -*-
"""

"""

import tifffile

import numpy as np

from scipy import ndimage
from skimage.morphology import disk
from skimage.segmentation import find_boundaries
from tqdm import tqdm

def gen_anno_bound(anno_im_f, dilation_iter=1):
    
    anno_im = tifffile.imread(anno_im_f)
    anno_im = np.moveaxis(anno_im, 0, 2)
    
    anno_bound_im = np.zeros(np.shape(anno_im), dtype=np.uint16)
    for i in tqdm(range(np.shape(anno_im)[2])):
        anno_bound_im[:,:,i] = find_boundaries(anno_im[:,:,i])
        
    for i in tqdm(range(np.shape(anno_bound_im)[2])):
        anno_bound_im[:,:,i] = ndimage.binary_dilation(anno_bound_im[:,:,i], 
                                                       structure=disk(1),
                                                       iterations=dilation_iter)
    
    return anno_bound_im

