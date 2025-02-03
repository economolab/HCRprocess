Formatting Inputs
------------------------------

Create an experiment directory in the following format:

YYYY-MM-DD_EXPERIMENT

Examples:

.. code-block:: 
   2024-12-20_FN-SNR-8
   2025-01-06_XII-2
   2023-06-08_MC-SC-5

Create a subdirectory within this folder for each round (pre-imaging round number is 0) in the following format:
```
YYYY-MM-DD_r#_HCR

Examples:
2024-12-23_r1_HCR
2025-01-06_r0_HCR
2023-06-14_r2_HCR
```
The ``HCR`` tag is required for HCRprocess to recognize the folder as a round folder. Within each round folder, place your nd2 files inside it in the following format:
```
UNIQUE-ID__gene1chan#_gene2chan#_gene3chan#_gene4chan##__HCR_EXTRA-INFO

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

