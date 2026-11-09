# -*- coding: utf-8 -*-
"""
Created on Wed Oct 28 12:56:04 2026

@author: jpv88
"""

# %% imports

import os
import tifffile

import numpy as np
import pandas as pd

from tqdm import tqdm

from analysis_utils import matrix_utils, filter_data, img_utils

#%% explicit and derived params

data_dir = 'D:\\2026-09-12_FN_SNr_2\\post\\core_output\\s01L'

uniq_id, post_dir, exp_name = matrix_utils.derive_params(data_dir)
data_f, reg_f, bin_f, mask_f, mask_im_f = matrix_utils.build_paths(data_dir, 
                                                                   uniq_id, 
                                                                   post_dir, 
                                                                   exp_name,
                                                                   mode='curated')

#%%

TF = [matrix_utils.check_for_tokens(f,['r1']) for f in reg_f]
r1_reg_im_f = os.path.join(data_dir, reg_f[TF][0])



