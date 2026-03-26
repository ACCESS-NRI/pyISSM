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

pyISSM conda environment
-------------------------------------------
An `environment.yml` file is provided in the repository to create an optional conda environment with all necessary dependencies to use pyISSM. You can create the environment using the following command:

.. code-block:: bash

   conda env create -f environment.yml

Once created, you can activate the environment with:

.. code-block:: bash

   conda activate pyissm

Alternatively, for NCI Gadi users, you can use the ACCESS-MED Containerised Conda Environment `analysis3` conda environment maintained within the `xp65` project. More information and instructions can be found on 
the `ACCESS-Analysis-Conda GitHub repository <https://github.com/ACCESS-NRI/ACCESS-Analysis-Conda>`_.

Verifying the Installation
--------------------------
After installation, you can verify the selected pyISSM version is installed correctly by running the following commands in a Python shell:

.. code-block:: python

   import pyissm
   print(pyissm.__version__)