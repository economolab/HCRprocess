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

#. Pitch and yaw of the blocking cut for each image file. These can be edited manually here and will automatically save. They are also loaded automatically into this table when saved within the blocking cut GUI. These values are auto loaded into the CCF registration GUI for each image file when it's opened. Pitch should be identical for each image file. The magnitude of the yaw should be identical as well, but its sign may need to be flipped depending on the orientation of the image file compared to whatever was used to calculate the pitch and yaw originally. This is discussed in more detail in the blocking cut GUI documentation. 

#. CCF registration functions. Press ``Launch block cut GUI`` to launch the blocking cut angle finding GUI. Files to use are selected within the GUI. This GUI only needs to be used once per experiment. Press ``CCF register this file`` to open the CCF registration GUI for the file selected in the table below. This GUI needs to be used once for every HCR image file in an experiment.

#. This table stores info about CCF registration completion for all image files. The ``CCF Registration Settings`` column shows whether a CCF registration settings file has been saved for this image file yet. This file contains information about fiducials, blocking cut angle, and 3D crop dimensions. The ``CCF Registered`` column shows whether the actual registration procedure has been successfully carried out using a settings file.

Cellpose-SAM Mask Segmentation and Quality Control
~~~~~~~~~~~~~~~~~~~~~~~

Final Cell Quality Control
~~~~~~~~~~~~~~~~~~~~~~~
