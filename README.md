# Installation

[Install MATLAB](https://www.mathworks.com/help/install/ug/install-products-with-internet-connection.html)

Clone this repository and run `HCRprocess.mlapp`

## ImageJ in MATLAB
Required for stitching, registration, and final processing

[Install Fiji](https://imagej.net/software/fiji/downloads)

[Add the ImageJ-MATLAB update site in ImageJ and increase your Java heap memory size in MATLAB](https://imagej.net/scripting/matlab)

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
```

# User Guide

# TODO

1) Add user guide
2) lipo gone
3) turn on and off blocks of processing pipeline
4) cellpose 3 model training and implementation
5) Batch registration
6) do all button
7) memory monitor
8) memory usage reduction and clearing
9) separate back end from gui for reg and unmix, run with progress dialog box
10) delete intermediates
11) delete or properly store files created by registration
