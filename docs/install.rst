Installation
------------------------------

`Install MATLAB <https://www.mathworks.com/help/install/ug/install-products-with-internet-connection.html>`_

Clone this repository and run ``HCRprocess.mlapp``

ImageJ in MATLAB: Required for stitching and final processing

`Install Fiji <https://imagej.net/software/fiji/downloads>`_

`Add the ImageJ-MATLAB update site in ImageJ and increase your Java heap memory size in MATLAB <https://imagej.net/scripting/matlab>`_

`How to increase Java heap memory size beyond the default maximum allowed value <https://www.mathworks.com/matlabcentral/answers/92813-how-do-i-increase-the-heap-space-for-the-java-vm-in-matlab#answer_183274>`_ A minimum of at least two to three times the size of the largest image you intend to work with is recommended. 

Set ``fijiSciptsDir`` in ``params.m`` to the correct Fiji scripts directory on your system. If you don't know where this is, you can use ``File -> Show Folder -> ImageJ`` in Fiji. 

Python in MATLAB: Required for registration and segmentation

`Install Anaconda Distribution <https://docs.anaconda.com/anaconda/install/>`_

Create an environment using Anaconda Navigator. The environment can be named whatever you want, but if you don't use the default name of ``HCRprocess``, be sure to change ``HCRprocessEnv`` in ``params.m`` to the correct name. 

In Anaconda Prompt, activate the environment, and then install the required packages:

.. code-block:: python

   conda activate HCRprocess
   pip install itk
   pip install itk-elastix
   pip install cellpose[gui]


Make sure you have `added your Anaconda scripts directory to your Path <https://www.architectryan.com/2018/03/17/add-to-the-path-on-windows-10/>`_, and that you have run ``conda init`` in your terminal at least once. ``HCRprocess`` calls Python using MATLAB's `system <https://www.mathworks.com/help/matlab/ref/system.html>`_ function and the `python <https://docs.python.org/3/using/cmdline.html>`_ command line invocation. 
