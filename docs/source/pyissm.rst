.. _api:

pyissm API Reference
===================
.. automodule:: pyissm
.. currentmodule:: pyissm

pyISSM comprises six major sub-packages:
- ` ``pyissm.model`` <model-class_>`_ - contains the central functionality, including the core ``Model`` class and all model `sub-classes <model-subclasses_>`_.
- ` ``pyissm.tools`` <tools_>`_ - contains diverse tools for interacting with the ISSM Model object.
- ` ``pyissm.data`` <data_>`_ - contains tailored tools for interacting with external datasets, including interpolation routines.
- ` ``pyissm.plot`` <plot_>`_ - contains tailored tools for visualising ISSM models.
- ` ``pyissm.analysis`` <analysis_>`_ - contains tailored post-processing tools.

.. _model-class:

pyISSM Model Class
-------------------
``pyISSM`` is developed around a central ``Model()`` class.

.. autosummary::
   :toctree: generated
   :recursive:

   pyissm.model.Model

.. _model-subclasses:

pyISSM Model sub-classes
-------------------
Each ISSM model contains a series of sub-classes that define model components and parameters (e.g., ``mesh``, ``geometry``, etc.).

.. autosummary::
   :toctree: generated
   :recursive:

   pyissm.model.classes

.. _tools:

Tools
-------------------
Utilities and helper functions for interacting with the ISSM ``Model`` object.

.. autosummary::
   :toctree: generated
   :recursive:

   pyissm.tools

.. _data:

Data
-------------------
Function for interacting with external datasets and interpolation routines.

.. autosummary::
   :toctree: generated
   :recursive:

   pyissm.data

.. _plot:

Plot
-------------------
Visualisation utilities, including plotting of model meshes, boundary conditions, fields, and results.

.. autosummary::
   :toctree: generated
   :recursive:

   pyissm.plot

.. _analysis:

Analysis
-------------------
Post-processing tools for analysing ISSM models, including derived fields and metrics.

.. autosummary::
   :toctree: generated
   :recursive:

   pyissm.analysis
