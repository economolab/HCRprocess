Input Files
------------------------------

Create an experiment directory in the following format:

YYYY-MM-DD_EXPERIMENT

Examples:

.. code-block::

   2024-12-20_FN-SNR-8
   2025-01-06_XII-2
   2023-06-08_MC-SC-5

Create a subdirectory within this folder for each round (pre-imaging round number is 0) in the following format:

YYYY-MM-DD_r#_HCR

Examples:

.. code-block:: 

   2024-12-23_r1_HCR
   2025-01-06_r0_HCR
   2023-06-14_r2_HCR

The ``HCR`` tag is required for HCRprocess to recognize the folder as a round folder. Within each round folder, place your nd2 files in the following format: 

UNIQUE-ID__gene1chan#_gene2chan#_gene3chan#_gene4chan#__HCR-EXTRA-INFO

Examples:

.. code-block::

   s04__eGFP488_dTom546_NT647__HCR-40x.nd2
   s03-left__Phox2b594_Snap25488_Calb1546_Slc32a1647__HCR-40x-zstack.nd2
   s01-ant__Phox2b594_Snap25488_Calb1546_Slc32a1647__HCR.nd2

Note the double underscores separating both the unique identifier from channels and channels from extra information. Make sure a single underscore separates each channel, and that channels are listed in the order they appear in the file (the order in which they are acquired in NIS-Elements). The unique identifier and extra information can be whatever you want, but don't include underscores in them. The ``HCR`` tag is required for HCRprocess to recognize a file as an HCR image file that requires processing. Here's an example of what a full experiment directory might look like, before any processing in ``HCRprocess``:

.. code-block::

   2024-11-10_MC-SC-4
   |  2024-11-10_r0_HCR
   |  |  s01-20x-overview.nd2*
   |  |  s02-20x-overview.nd2*
   |  |  s01-ROI.png*
   |  |  s02-ROI.png*
   |  |  s01__eGFP488_dTom546_NT647__HCR.nd2
   |  |  s02__eGFP488_dTom546_NT647__HCR.nd2
   |  2024-11-13_r1_HCR
   |  |  s01-20x-overview.nd2*
   |  |  s02-20x-overview.nd2*
   |  |  s01__Phox2b594_Snap25488_Calb1546_Slc32a1647__HCR-zstack.nd2
   |  |  s02__Phox2b594_Snap25488_Calb1546_Slc32a1647__HCR-zstack.nd2
   |  2024-11-16_r2_HCR
   |  |  s01-20x-overview.nd2*
   |  |  s02-20x-overview.nd2*
   |  |  s01__Nr4a2594_Snap25488_Zfhx4546_Ebf3647__HCR-40x.nd2
   |  |  s02__Nr4a2594_Snap25488_Zfhx4546_Ebf3647__HCR-40x.nd2

* = ignored by automatic file detection

``HCRprocess`` can also autodetect and process histology image files. First add a histology subdirectory in the experiment folder using the following formatting:

YYYY-MM-DD_histo

The ``histo`` tag is required for HCRprocess to recognize the folder as a histology folder. Within the histology folder, place your image files in the following format:

REGION#__fluorophore1chan#_fluorophore2chan#_fluorophore3chan#__histo-extra-info

Examples: 

.. code-block::

   MC2__dTomTRITC_NTCy5__histo.vsi
   SC5__eGFPFITC_dTomTRITC_NTCy5__histo-10x.vsi
   FN1__eGFPFITC_NTCy5__histo.nd2

Make sure a single underscore separates each channel, and that channels are listed in the order they appear in the file (the order in which they are acquired). The ``histo`` tag is required for HCRprocess to recognize a file as a histology image file that requires processing. Currently ``HCRprocess`` can automatically recognize and process histology image files of the following formats: vsi, nd2. 

