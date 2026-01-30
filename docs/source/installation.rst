Installation
============

.. warning::
   This is an experimental release. Documentation remains a work in progress. Some sections may be incomplete or under development.

pyISSM can be installed using `pip` or `conda`, or by cloning the repository from GitHub. Below are instructions for each method.

Installing via pip
-----------------
You can install pyISSM using `pip`, the Python package manager.
See the official release on PyPI: https://pypi.org/project/pyissm/

Run the following command in your terminal:

.. code-block:: bash

   pip install pyissm

Installing via conda
-------------------
If you prefer using `conda`, you can install pyISSM from the `accessnri` channel.  
https://anaconda.org/channels/accessnri/packages/pyissm/overview

.. code-block:: bash

   conda install accessnri::pyissm

Installing from GitHub (Development Version)
-------------------------------------------
You can also install the latest development version directly from the GitHub repository:

1. Clone the repository:

   .. code-block:: bash

      git clone https://github.com/ACCESS-NRI/pyISSM.git

2. Navigate to the cloned directory:

   .. code-block:: bash

      cd pyissm

3. Install pyISSM using pip:

   .. code-block:: bash

      pip install .

.. note::
   Installing from GitHub is recommended only if you want the latest features or are contributing to pyISSM development.

Verifying the Installation
--------------------------
After installation, you can verify that pyISSM is installed correctly by running the following commands in a Python shell:

.. code-block:: python

   import pyissm
   print(pyissm.__version__)