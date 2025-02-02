Manual GUIs
------------------------------

Multiple operations in HCRprocess require manual intervention. These operations call popup external GUIs, and instructions for using these GUIs are described below. 

Spectral Unmixing GUI
~~~~~~~~~~~~~~~~~~~~~~~

This GUI is called as part of the spectral unmixing operation. It is a GUI developed by Michael Economo for the purpose of manually setting parameters for linear unmixing of a multispectral image. 

Basic usage instructions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. image:: doc_assets/unmixing_manual_gui_final.png
    :width: 800px
    :align: center
    :alt: unmixing manual GUI

#. Displays the current image being unmixed as well as its name. Shows only one slice from the middle of the z-stack. 

#. Channels panel: displays all channels contained in the image. Turn channels on and off using the `On` checkbox. Use the `Auto` radio button to enable auto contrast scaling, or disable it to manually set pixel minimum and maximum values. Change the color of a channel by clicking the color indicator. Check the `HiLo LUT` checkbox to enable HiLo pixel coloring. Channels are displayed here from left to right and can be indexed from 1 to the total number of channels. Those indices are the ones used in the spectral unmixing table.

#. Spectral unmixing table: This table is where spectral unmixing parameters are manually set. `Unmix from` denotes the channel whose signal will be subtracted from the target `Unmix to` channel. For instance if channel 2 contained spectral bleed from fluorophores in channel 3 that you wanted to remove, then channel 3 would be your `Unmix from` channel, and channel 2 would be your `Unmix to` channel. `Scale` and `Offset` are the input linear weights for an unmixing operation. `Enable` determines whether that particular linear unmixing will be carried out. Unmixing operations are carried out in the order they are listed in the table, from top to bottom. 

#. Outdated lipofuscin removal model that will eventually be removed

#. Image shifting panel: can shift individual channels a desired number of pixels in-plane. This is helpful for imaging systems where certain channels are offset a pixel or two from others, however, this is not really an issue on the Nikon spinning disk confocal. 

#. Press `Update Image` to apply transformations set in other GUI components. Press `Save settings` when you are satisfied with the unmixing parameters that have been set. This will also close the GUI. 

Advanced usage instructions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. image:: doc_assets/unmixing_examples_final.png
    :width: 800px
    :align: center
    :alt: unmixing manual GUI


Finding Fiducials GUI
~~~~~~~~~~~~~~~~~~~~~~~

Lipofuscin Removal GUI
~~~~~~~~~~~~~~~~~~~~~~~

