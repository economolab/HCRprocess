HCRprocess is an HCR-mFISH image processing pipeline written in MATLAB and Python 3 by Jack Vincent. It is composed of components developed by Jack Vincent, Michael Economo, and Will Cunningham.   

![HCR example](HCRprocess/doc_assets/beauty_HCR.png)

# TODO

1) Add user guide
2) lipo gone
3) turn on and off blocks of processing pipeline
4) cellpose 3 model training and implementation
5) Batch registration
6) do all button
8) memory usage reduction and clearing
9) separate back end from gui for reg and unmix, run with progress dialog box
10) delete intermediates
11) delete or properly store files created by registration
12) fix tiff writing in MATLAB (and tiff reading, need dedicated functions for both that are called for everything, currently using hacky solution of reopening and writing in imagej for all writes, we already have write while incrementing ifd functionality, just need read while incrementing ifd functionality)
