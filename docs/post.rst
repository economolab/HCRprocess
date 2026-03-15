Post-processing
------------------------------

The final output of the core processing pipeline for each image file is a folder containing every channel acquired across all rounds as its own separate tif, with all images spectrally unmixed and registered to a single registration round. This is stored in the core_output folder.

.. code-block:: 

   experiment_directory
   |  post
   |  |  core_output

The purpose of the post-processing pipeline is to convert this folder into a rigorously quality controlled cells x genes/features matrix for subsequent data analysis.

Allen CCF Registration
~~~~~~~~~~~~~~~~~~~~~~~

This tab is where image files are registered using an affine fiducial-based warp to the Allen Mouse Brain Common Coordinate Framework. 

.. image:: doc_assets/ccf_gui_final.png
    :width: 800px
    :align: center
    :alt: unmixing GUI

1. Pitch and yaw of the blocking cut for each image file. These can be edited manually here and will automatically save. They are also loaded automatically into this table when saved within the blocking cut GUI. These values are auto loaded into the CCF registration GUI for each image file when it's opened. Pitch should be identical for each image file. The magnitude of the yaw should be identical as well, but its sign may need to be flipped depending on the orientation of the image file compared to whatever was used to calculate the pitch and yaw originally. This is discussed in more detail in the blocking cut GUI documentation. 

Cellpose-SAM Mask Segmentation and Quality Control
~~~~~~~~~~~~~~~~~~~~~~~

Final Cell Quality Control
~~~~~~~~~~~~~~~~~~~~~~~
