.. _api:

pyissm API Reference
===================
.. automodule:: pyissm
.. currentmodule:: pyissm

pyISSM comprises five major sub-packages:

- ``pyissm.model`` - contains the central functionality, including the core ``Model`` class and all model sub-classes.
- ``pyissm.tools`` - contains diverse tools for interacting with the ISSM Model object.
- ``pyissm.data`` - contains tailored tools for interacting with external datasets, including interpolation routines.
- ``pyissm.plot`` - contains tailored tools for visualising ISSM models.
- ``pyissm.analysis`` - contains tailored post-processing tools.

.. _model-class:

Model Class
-------------------
``pyISSM`` is developed around a central ``Model()`` class.

.. autosummary::
   :toctree: api
   :recursive:

   pyissm.model.Model

.. _model-subclasses:

Model sub-classes
-------------------
Each ISSM model contains a series of sub-classes that define model components and parameters (e.g., ``mesh``, ``geometry``, etc.).

.. autosummary::
   :toctree: api
   :recursive:

   pyissm.model.classes

.. _model-operations:

Model operations
-------------------
Each ISSM model is constructed and parameterised using various model operations.

.. autosummary::
   :toctree: api
   :recursive:

   pyissm.model.bc
   pyissm.model.execute
   pyissm.model.io
   pyissm.model.mesh
   pyissm.model.param

.. _tools:

Tools
-------------------
Utilities and helper functions for interacting with the ISSM ``Model`` object.

.. autosummary::
   :toctree: api
   :recursive:

   pyissm.tools

.. _data:

Data
-------------------
Function for interacting with external datasets and interpolation routines.

.. autosummary::
   :toctree: api
   :recursive:

   pyissm.data

.. _plot:

Plot
-------------------
Visualisation utilities, including plotting of model meshes, boundary conditions, fields, and results.

.. autosummary::
   :toctree: api
   :recursive:

   pyissm.plot

.. _analysis:

Analysis
-------------------
Post-processing tools for analysing ISSM models, including derived fields and metrics.

.. autosummary::
   :toctree: api
   :recursive:

   pyissm.analysis
