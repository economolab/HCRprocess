Manual GUIs
------------------------------

Multiple operations in HCRprocess require manual intervention. These operations call popup external GUIs, and instructions for using these GUIs are described below. 

Spectral Unmixing GUI
~~~~~~~~~~~~~~~~~~~~~~~

This GUI is called as part of the spectral unmixing operation. It is a GUI developed by Michael Economo for the purpose of manually setting parameters for linear unmixing of a multispectral image. See `this article from Nikon <https://www.microscopyu.com/techniques/confocal/spectral-imaging-and-linear-unmixing>`_ for a general overview of spectral imaging and linear unmixing. 

Basic usage instructions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. image:: doc_assets/unmixing_manual_gui_final.png
    :width: 800px
    :align: center
    :alt: unmixing manual GUI

#. Displays the current image being unmixed as well as its name. Shows only one slice from the middle of the z-stack. 

#. Channels panel: displays all channels contained in the image. Turn channels on and off using the ``On`` checkbox. Use the ``Auto`` radio button to enable auto contrast scaling, or disable it to manually set pixel minimum and maximum values. Change the color of a channel by clicking the color indicator. Check the ``HiLo LUT`` checkbox to enable HiLo pixel coloring. Channels are displayed here from left to right and can be indexed from 1 to the total number of channels. These indices are the ones used in the spectral unmixing table.

#. Spectral unmixing table: This table is where spectral unmixing parameters are manually set. ``Unmix from`` denotes the channel whose signal will be subtracted from the target ``Unmix to`` channel. For instance if channel 2 contained spectral bleed from fluorophores in channel 3 that you wanted to remove, then channel 3 would be your ``Unmix from`` channel, and channel 2 would be your ``Unmix to`` channel. ``Scale`` and ``Offset`` are the input linear weights for an unmixing operation. ``Enable`` determines whether that particular linear unmixing will be carried out. Unmixing operations are carried out in the order they are listed in the table, from top to bottom. 

#. Outdated lipofuscin removal model that will eventually be removed

#. Image shifting panel: can shift individual channels a desired number of pixels in-plane. This is helpful for imaging systems where certain channels are offset a pixel or two from others, however, this is not really an issue on the Nikon spinning disk confocal. 

#. Press ``Update Image`` to apply transformations set in other GUI components. Press ``Save settings`` when you are satisfied with the unmixing parameters that have been set. This will also close the GUI. 

Advanced usage instructions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The goal of spectral unmixing is to remove all spectral bleed between channels. There are only two failure modes: under-unmixing and over-unmixing. Below is a detailed description of how to unmix an HCR-mFISH multispectral image while avoiding these two scenarios.

.. image:: doc_assets/unmixing_examples_final.png
    :width: 800px
    :align: center
    :alt: unmixing manual GUI

* Step 1 is determining which channels need to be unmixed from which other channels, as well as the order in which these operations will be carried out. Generally fluorophores will only bleed into spectrally adjacent channels, e.g. 488 might bleed into 514, but not 546, 594, or 647. It's typically a good idea to just unmix all spectrally adjacent channels from one another even if the magnitude of some of those unmixing operations is miniscule. Two considerations need to be taken into account when determining unmixing order: it's a good idea to unmix channels with strong signal before channels with weak signal, and it's a good idea to unmix channels with only one possible contamination source before channels with two possible contamination sources. The reasoning behind these principles is similar, unmixing is a destructive operation that will inevitably remove some true signal in addition to spectral bleed. Channels that are unmixed first will have more true signal removed, since their signal hasn't yet been removed from the spectral bleed source channel being subtracted from them. As an example, when unmixing 546, 594, and 647: only 594 needs to be unmixed twice, and channel signal strength is rougly 647 > 594 > 546. In this case the most logical order would be 594 from 647, 594 from 546, 546 from 594, and lastly 647 from 594. 
* Step 2 is actually setting the parameters for unmixing one channel from another. In general, an offset is not needed, except in cases of extreme spectral bleed. A good starting point is to identify neurons that contain solely or near solely shadow bleed from adjacent channels. Ideally the channel should be unmixed to the point that these shadow neurons are no longer visible, and just look like background. An easy way to find neurons like this is to intially heavily oversubtract and then look for neuron shadows full of desaturated pixels. Then, initially coarsely and then later gently, increment the unmixing scale downward until desaturated pixels disappear and the neuron becomes indistinguishable from background (this means its neither noticeably darker nor noticeably brighter than tissue background pixels). An initial increment of 0.1-0.2 down to a final fine-tuning increment of 0.01-0.05 works well. Check many areas of the slice to ensure the parameters you have set look correct everywhere. It can often be helpful to flick the spectral bleed source channel on and off to see what it is you need to unmix from the target channel, as well as determine which neurons are genuine signal and which are spectral bleed shadow. 


Finding Fiducials GUI
~~~~~~~~~~~~~~~~~~~~~~~

This GUI is called as part of the registration operation. It is a GUI developed by Michael Economo for the purpose of manually marking fiducials for registration of two images. See `this Wikipedia article <https://en.wikipedia.org/wiki/Image_registration>`_ for a general overview of image registration, and `this webpage from the University of Edinburgh <https://homepages.inf.ed.ac.uk/rbf/HIPR2/affine.htm>`_ for a general overview of affine transformations. 

Basic usage instructions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. image:: doc_assets/find_fiducials_gui_final.png
    :width: 800px
    :align: center
    :alt: find fiducials manual GUI

#. Fixed image display: navigate through Z using the frame scrollbar, press the ``Auto`` radio button to enable automatic contrast adjustment, otherwise manually set pixel minimum and maximum values. The green dot is a movable entity that denotes where the next fiducial will be added. Red dots indicate previously marked fiducials.

#. Moving image display: usage similar as fixed image display. Movable fiducial dot is magenta. 

#. Axes controls: check ``Scroll both`` to scroll both axes in Z simultaneously (this checkbox is not currently functional and will just set the moving image Z-frame to the same index as the fixed image Z-frame). Check ``Overlay`` to directly overlay the fixed image in green on the moving image in magenta. Only usable if a transform has been performed. Also enables the green fixed image fiducial dot in the overlay image. All other controls in this panel are deprecated. 

#. Lists all marked fiducials. Click a row to select it. 

#. Control panel: press ``Add point`` to add the current locations of the fixed image green dot and moving image magenta dot as fiducials. Press ``Affine Warp`` to perform an affine warp using the current fiducials. This is useful for previewing current progress. ``b-spline Warp`` is currently deprecated. Press ``Delete point`` to delete the currently selected row in the fiducials table. Press ``Save points`` when you are satisfied with the fiducials that have been marked. This will also close the GUI. 

Advanced usage instructions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The goal of finding fiducials is to produce a high quality initial affine transform. HCRprocess's registration is two-step: with an initial fiducial-based affine transform followed by an intensity-based elastix B-spline warp. The intensity-based B-spline warp is very good at aligning individual cells and anatomical features, but it requires that the two input image volumes already be at least approximately aligned. Without a sufficient number of high quality fiducials, the initial affine transform will not be effective enough to align the slices, and the intensity-based B-spline warp will fail. Below are a number of tips to ensure registration success.

.. image:: doc_assets/find_fiducials_example.png
    :width: 800px
    :align: center
    :alt: unmixing manual GUI

* Example fiducial: these two image volumes were acquired using different fluorophores, one with NeuroTrace and the other with Snap25 fluorescence. Note the presence of multiple landmarks for alignment, including the cluster of 4 neurons on the right, as well as the elongated neuron on the left. Often the shape of a single neuron is not adequate to ensure fiducial alignment, make sure the area surrounding the fiducial lines up as well. Also note the presence of nuclei and non-neuronal labeling in the NeuroTrace sample. This is normal, but it does mean that certain features will be present in one stack but not the other, and it's up to the user to identify which features are missing due to different labeling strategies and which are absent due to fiducial misalignment.  
* Always ensure adequate spread of fiducials. The goal of the initial affine transform is not to align one area of the slice very well, but rather to align the entire slice with moderate accuracy. Do not concentrate fiducials in-plane, make sure you are evenly sampling the dorsal, ventral, medial, and lateral portions of the slice. Also be sure to put more fiducials on the z-extremes (superficial and deep) of the slice. These areas are the most prone to misalignment, and aligning both helps ensure the proper alignment of the intervening tissue. In general, 1/3 of fiducials should be concentrated in the superficial edge of the slice, 1/3 in the deepest few z-planes of the slice, and 1/3 in all the intervening z-planes. 
* Aligning images from pre-HCR imaging to post-HCR imaging is much more difficult than aligning post to post, since the initial HCR treatments significantly shrink and deform the tissue section, and pre-HCR treatments are captured with NeuroTrace as opposed to Snap25 fluorescence. It is recommended to mark at least 15-20 fiducials for these transforms and perhaps up to as much as 25-30 for highly deformed sections. For post-to-post registration, 10-15 fiducials is generally adequate, although as many as 20 may be needed if imaging alignment is poor. 
* The affine transform is fundamentally limited in its ability to capture non-linear deformations. If adding more fiducials is not improving image alignment, this may simply be a limitation of the affine transform, and further alignment must be delegated to the intensity-based B-spline warp. 

Lipofuscin Removal GUI
~~~~~~~~~~~~~~~~~~~~~~~

