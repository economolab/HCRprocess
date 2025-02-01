Main
------------------------------

The ``Main`` tab is the central tab of HCRprocess from which all files in a round can be reviewed. Beneath the main tab's file tables are main functions which can be used in any tab of HCRprocess.

.. image:: doc_assets/main_GUI_final.png
    :width: 1200px
    :align: center
    :alt: main GUI

#. Click on each tab to navigate around HCRprocess
#. Core processing pipeline table: gives the status of each core processing operation for each file in the currently selected round. ✅ Green checkmarks denote completed operations, ❌ red x's denote incomplete operations, and 🗿 gray Moais denote operations that have been disabled for this round. 
#. Manual settings table: gives the status of each manual settings operation for each file in the currently selected round.
#. Post processing table: gives the status of each post processing operation for this experiment. Post processing operations are not round-specific; they are completed only once per experiment. 
#. Experiment directory: Shows the current experiment directory as well as each round folder contained within the directory. Click on round folders to change the active round folder. Click ``Change Experiment Directory`` to change the experiment directory using a folder selection dialog box. 
#. Processing pipelines: Switches which can enable and disable portions of the core and post-processing pipelines. The ``Core Processing Pipeline`` switches are round-specific. ``Stitching`` is a required operation before all subsequent operations, therefore it is permanently enabled for all rounds.
#. Main functions:
   * Registration Round: Set the registration round using a drop-down menu. This is the round that all other rounds will be registered to. Its registration operation is automatically disabled, since there's no need to             register a round to itself. Upon loading an experiment directory for the first time, HCRprocess will attempt to find the registration round by looking for a registration round token ('r1' by default). If it's unable to       find it, this parameter must be set manually before carrying out registration operations. 
