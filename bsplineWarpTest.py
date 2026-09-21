# -*- coding: utf-8 -*-
"""
Created on Mon Aug 24 17:19:32 2026

@author: jpv88
"""

from scipy.io import loadmat
from scipy.ndimage import affine_transform

import numpy as np
import SimpleITK as sitk

import tifffile

#%%

f_reg_data = 'D:\\2026-01-16_MC_SC_17\\post\\CCF_vols\\s03L\\s03LCCFReg.mat'
f_cropped_ccf_vol = 'D:\\2026-01-16_MC_SC_17\\post\\CCF_vols\\s03L\\s03L_cropped_ccf_vol.mat'
f_affine_ccf_vol = 'D:\\2026-01-16_MC_SC_17\\post\\CCF_vols\\s03L\\s03Lannotation.tif'

#%%

data = loadmat(f_reg_data, struct_as_record=False, squeeze_me=True)
appStateStruct = data['appStateStruct']

pointsData = appStateStruct.pointsData
pointsCCF = pointsData.ccf
pointsCCFWarp = pointsData.warp_ccf
pointsSamp = pointsData.samp
sampDirection = appStateStruct.sampDirection
sampOrigin = appStateStruct.sampOrigin
sampSize = appStateStruct.sampSize
sampSpacing = appStateStruct.sampSpacing

data = loadmat(f_cropped_ccf_vol, struct_as_record=False, squeeze_me=True)
ccf_vol = data['anno']

#%%

ccf_vol_af = tifffile.imread(f_affine_ccf_vol)
ccf_vol_af = np.moveaxis(ccf_vol_af, 0, 2)

#%%

fixed_reference = sitk.Image([int(s) for s in sampSize], sitk.sitkUInt8)
fixed_reference.SetSpacing([float(s) for s in sampSpacing])
fixed_reference.SetOrigin([float(o) for o in sampOrigin])
fixed_reference.SetDirection([float(d) for d in sampDirection.flatten()])

mesh_size = [1, 1, 1]
bspline_transform = sitk.BSplineTransformInitializer(fixed_reference, mesh_size)

landmark_initializer = sitk.LandmarkBasedTransformInitializerFilter()

fixed_points = [float(v) for v in pointsSamp.flatten()]
# moving_points = [float(v) for v in pointsCCF.flatten()]
moving_points = [float(v) for v in pointsCCFWarp.flatten()]
landmark_initializer.SetFixedLandmarks(fixed_points)
landmark_initializer.SetMovingLandmarks(moving_points)

landmark_initializer.SetReferenceImage(fixed_reference)
bspline_transform = landmark_initializer.Execute(bspline_transform)

#%%

from scipy.interpolate import RBFInterpolator
from scipy.ndimage import map_coordinates


tps = RBFInterpolator(pointsSamp, pointsCCFWarp, kernel='thin_plate_spline', smoothing=50)

out_shape = sampSize

# --- out_shape stays in (X, Y, Z) convention ---
out_shape = tuple(int(s) for s in sampSize)  # (X, Y, Z) = (5170, 5170, 40)
Nx, Ny, Nz = out_shape

# --- Output array, same (X, Y, Z) convention as ccf_vol ---
registered = np.zeros((Nx, Ny, Nz), dtype=ccf_vol_af.dtype)

for z in range(Nz):
    xx, yy = np.meshgrid(np.arange(Nx), np.arange(Ny), indexing='ij')
    zz = np.full_like(xx, z)
    query_pts = np.stack([xx.ravel(), yy.ravel(), zz.ravel()], axis=-1)  # (x, y, z)

    sampled_xyz = tps(query_pts)  # (M, 3), (x, y, z) order, in ccf_vol space

    coords = np.stack([
        sampled_xyz[:, 0].reshape(Nx, Ny),  # x -> ccf_vol axis 0
        sampled_xyz[:, 1].reshape(Nx, Ny),  # y -> ccf_vol axis 1
        sampled_xyz[:, 2].reshape(Nx, Ny),  # z -> ccf_vol axis 2
    ], axis=0)

    registered[:, :, z] = map_coordinates(ccf_vol_af, coords, order=1, mode='constant', cval=0.0)

    print(f"slice {z}/{Nz}")
    
#%%

to_save = np.moveaxis(registered, 2, 0)
tifffile.imwrite('test.tiff', to_save)


#%%

af = appStateStruct.af  # 4x4, row-vector convention

M = af[:3, :3].T   # 3x3 linear part, transposed for column-vector use
t = af[3, :3]      # translation vector

# Now: p_out (column vector) = M @ p_in + t

M_inv = np.linalg.inv(M)
t_inv = -M_inv @ t

out_shape_xyz = tuple(int(s) for s in sampSize)  # (X, Y, Z), fixed/Samp dims

registered = affine_transform(
    ccf_vol,
    matrix=M_inv,
    offset=t_inv,
    output_shape=out_shape_xyz,
    order=1,
    mode='constant',
    cval=0.0
)