# Installation

[Install MATLAB](https://www.mathworks.com/help/install/ug/install-products-with-internet-connection.html)

Clone this repository and run `HCRprocess.mlapp`

## ImageJ in MATLAB
Required for stitching and final processing

[Install Fiji](https://imagej.net/software/fiji/downloads)

[Add the ImageJ-MATLAB update site in ImageJ and increase your Java heap memory size in MATLAB](https://imagej.net/scripting/matlab)

[How to increase Java heap memory size beyond the default maximum allowed value](https://www.mathworks.com/matlabcentral/answers/92813-how-do-i-increase-the-heap-space-for-the-java-vm-in-matlab#answer_183274)

Set `fijiSciptsDir` in `params.m` to the correct Fiji scripts directory on your system. If you don't know where this is, you can use `File -> Show Folder -> ImageJ` in Fiji. 

## Python in MATLAB
Required for registration and segmentation

[Install Anaconda Distribution](https://docs.anaconda.com/anaconda/install/)

Create an environment using Anaconda Navigator. The environment can be named whatever you want, but if you don't use the default name of `HCRprocess`, be sure to change `HCRprocessEnv` in `params.m` to the correct name. 

In Anaconda Prompt, activate the environment, and then install the required packages:

```
conda activate HCRprocess
pip install itk
pip install itk-elastix
pip install cellpose[gui]
```

Set `anacondaSciptsDir` in `params.m` to the correct Anaconda scripts directory on your system. `HCRprocess` calls Python using MATLAB's [system](https://www.mathworks.com/help/matlab/ref/system.html) function and the [python](https://docs.python.org/3/using/cmdline.html) command line invocation. 

# User Guide

## Input File Structuring
Create an experiment directory in the following format:
```
YYYY-MM-DD_EXPERIMENT

Examples:
2024-12-20_FN-SNR-8
2025-01-06_XII-2
2023-06-08_MC-SC-5
```
Create a subdirectory within this folder for each round (pre-imaging round number is 0) in the following format:
```
YYYY-MM-DD_r#

Examples:
2024-12-23_r1
2025-01-06_r0
2023-06-14_r2
```
Within each round folder, create a subdirectory called `raw` and place your nd2 files inside it in the following format:
```
UNIQUE-ID__gene1chan#_gene2chan#_gene3chan#_gene4chan##__EXTRA-INFO

Examples:
s04__eGFP488_dTom546_NT647__40x-zstack.nd2
s03-left__Phox2b594_Snap25488_Calb1546_Slc32a1647__40x-zstack.nd2
s01-ant__Phox2b594_Snap25488_Calb1546_Slc32a1647__20x-reg.nd2
```
Note the double underscores separating both the unique identifier from channels and channels from extra information. Make sure a single underscore separates each channel, and that channels are listed in the order they appear in the file (the order in which they are acquired in NIS-Elements). The unique identifier and extra information can be whatever you want, but don't include underscores in them. Here's an example of what a full experiment directory might look like, before any processing in `HCRprocess`:
```
2024-11-10_MC-SC-4
|  2024-11-10_r0
|  | raw
|  |  |  s01-20x-overview.nd2*
|  |  |  s02-20x-overview.nd2*
|  |  |  s01-ROI.png*
|  |  |  s02-ROI.png*
|  |  |  s01__eGFP488_dTom546_NT647__40x-zstack.nd2
|  |  |  s02__eGFP488_dTom546_NT647__40x-zstack.nd2
|  2024-11-13_r1
|  |  raw
|  |  |  s01-20x-overview.nd2*
|  |  |  s02-20x-overview.nd2*
|  |  |  s01__Phox2b594_Snap25488_Calb1546_Slc32a1647__40x-zstack.nd2
|  |  |  s02__Phox2b594_Snap25488_Calb1546_Slc32a1647__40x-zstack.nd2
|  2024-11-16_r2
|  |  raw
|  |  |  s01-20x-overview.nd2*
|  |  |  s02-20x-overview.nd2*
|  |  |  s01__Nr4a2594_Snap25488_Zfhx4546_Ebf3647__40x-zstack.nd2
|  |  |  s02__Nr4a2594_Snap25488_Zfhx4546_Ebf3647__40x-zstack.nd2

* = ignored by automatic file detection
```
Files containing any of the following case-insensitive keywords are ignored by the automatic file detection system: `ROI`, `overview`, `TileConfiguration`. Several subdirectories will be automatically created the first time `HCRprocess` is run with an experiment directory.   

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
