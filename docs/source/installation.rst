Installation
============

.. warning::
   This is an experimental release. Documentation remains a work in progress. Some sections may be incomplete or under development.


pyISSM can be installed using `pip` or `conda`, or by cloning the repository from GitHub. Below are instructions for each method.

## Installing via pip
To install pyISSM using pip, run the following command in your terminal:
   .. code-block:: bash
      pip install pyissm

## Installing via conda
If you prefer using conda, you can install pyISSM from the accessnri channel:
   .. code-block:: bash
      conda install accessnri::access-cryosphere-data-pool

## Installing from GitHub
The latest development version of pyISSM can be installed by cloning the GitHub repository.

1. Clone the repository:
   .. code-block:: bash
      git clone https://github.com/ACCESS-NRI/pyISSM.git

2. Navigate to the cloned directory:
   .. code-block:: bash
      cd pyissm
3. Install pyISSM using pip:
   .. code-block:: bash
      pip install .

## Verifying the Installation
After installation, you can verify that pyISSM is installed correctly by running the following command in a Python shell:
   .. code-block:: python
      import pyissm
      print(pyissm.__version__)