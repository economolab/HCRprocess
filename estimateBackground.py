# -*- coding: utf-8 -*-
"""
Use photutils to estimate the background of a 2D grayscale or 3D multichannel image.
Background is estimated using the sigma-clipped median in each box of a grid 
that covers the input data to create a low-resolution background map.
"""

import numpy as np

from astropy.stats import SigmaClip
from photutils.background import Background2D, MedianBackground


def estimateBackground(data):

    if data.ndim < 2 or data.ndim > 3:
        raise ValueError(
            "Input image must be either 2D for a single channel, or 3D for multichannel")

    if data.ndim == 2:
        data = np.expand_dims(data, axis=2)

    num_chans = np.shape(data)[2]
    
    box_size = 50
    filter_size = 3

    sigma_clip = SigmaClip(sigma=3.0)
    bkg_estimator = MedianBackground()

    background = np.zeros(np.shape(data))
    
    for i in range(num_chans):
        
        bkg = Background2D(data[:,:,i],
                           box_size,
                           filter_size=filter_size,
                           sigma_clip=sigma_clip,
                           bkg_estimator=bkg_estimator)
        
        background[:,:,i] = bkg.background
        
    return background

data = np.asarray(data)
background = estimateBackground(data)
