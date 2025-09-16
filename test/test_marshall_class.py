"""
Test marshall_class() functions on all classes within pyissm.param.XXX
"""

import inspect
import importlib
import pkgutil
import pytest
import numpy as np
import os
import tempfile

## Import pyISSM
import sys
sys.path.append('/Users/lawrence.bird/pyISSM/src/')
import pyissm

def build_base_model():
    """
    Create a base md model with default sub-classes.

    Returns
    -------
    md : pyissm.Model
        A configured model instance with mesh, friction, and geometry data.
    """

    ## Initialise model
    md = pyissm.Model()

    ## Define mesh information
    domain_name = os.path.join(os.path.dirname(__file__), 'assets/SquareIceShelf/DomainOutline.exp')
    rift_name = ''
    resolution = 100000
    area = resolution ** 2

    ## Create the mesh
    elements, x, y, segments, segmentmarkers = pyissm.utils.wrappers.Triangle(domain_name, rift_name, area)

    ## Assign mesh to model
    md.mesh.x = x
    md.mesh.y = y
    md.mesh.elements = elements
    md.mesh.segments = segments.astype(int)
    md.mesh.segmentmarkers = segmentmarkers.astype(int)
    md.mesh.numberofvertices = np.size(md.mesh.x)
    md.mesh.numberofelements = np.size(md.mesh.elements, axis = 0)
    md.mesh.vertexonboundary = np.zeros(md.mesh.numberofvertices, dtype = int)
    md.mesh.vertexonboundary[segments[:, 0:2] - 1] = 1
    md.mesh.vertexconnectivity = pyissm.utils.wrappers.NodeConnectivity(md.mesh.elements, md.mesh.numberofvertices)
    md.mesh.elementconnectivity = pyissm.utils.wrappers.ElementConnectivity(md.mesh.elements, md.mesh.vertexconnectivity)

    ## Assign defaults
    md.friction.coefficient = 100 * np.ones((md.mesh.numberofvertices))
    md.geometry.thickness = 100 * np.ones((md.mesh.numberofvertices))

    return md


def get_all_param_classes():
    """
    Discover and return all parameter classes from pyissm.param submodules.

    This function iterates through all submodules in the pyissm.param package
    and extracts classes that are defined within each module (not imported).

    Returns
    -------
    list of tuple
        A list of tuples containing (module_name, class_name, class_object)
        for each discovered parameter class.
    """

    classes = []
    package = pyissm.param
    for _, mod_name, ispkg in pkgutil.iter_modules(package.__path__):
        if ispkg:
            continue
        mod = importlib.import_module(f"{package.__name__}.{mod_name}")
        for name, cls in inspect.getmembers(mod, inspect.isclass):
            if cls.__module__ == mod.__name__:
                classes.append((mod.__name__, name, cls))
    return classes


## Run unit tests
all_results = []

@pytest.mark.parametrize("module_name,class_name,cls",
                         get_all_param_classes())
def test_marshall_class(module_name, class_name, cls):
    """
    Test marshall_class method for a parameter class.

    This test function creates a base model, instantiates the given parameter
    class, and tests its marshall_class method by attempting to write the
    class data to a temporary file.

    Parameters
    ----------
    module_name : str
        Full module name of the parameter class (e.g., 'pyissm.param.friction').
    class_name : str
        Name of the parameter class to test.
    cls : type
        The parameter class object to instantiate and test.

    Raises
    ------
    Exception
        If the marshall_class method fails during execution.

    Notes
    -----
    The test will skip classes that don't have a marshall_class method.
    For classes that fail the test, the exception is re-raised to ensure
    pytest properly reports the failure.
    """

    ## Create base model
    md = build_base_model()
    subclass_name = module_name.split('.')[-1]

    ## Attempt to get existing sub-class instance
    if hasattr(md, subclass_name):
        base_instance = getattr(md, subclass_name)
    else:
        ## Try to create a default instance if the attribute doesn't exist
        try:
            instance = cls()
            setattr(md, subclass_name, instance)
            base_instance = instance
        except (ModuleNotFoundError, AttributeError):
            ## Fallback if no default is available
            base_instance = None

    ## Create the class instance, inheriting from base_instance if possible
    try:
        if base_instance is not None:
            instance = cls(base_instance)
        else:
            instance = cls()
    except TypeError:
        ## Fallback if constructor does not accept a base instance
        instance = cls()

    ## Attach the class to the model
    setattr(md, subclass_name, instance)

    ## Print
    class_id = f"{module_name}.{class_name}"
    print(f"\nTesting class: {class_id}")

    ## If marshal_class() exists, test it.
    if hasattr(instance, "marshall_class"):
        try:
            prefix = "test_"
            with tempfile.NamedTemporaryFile() as fid:
                instance.marshall_class(fid, prefix, md)
            print(f"--> SUCCESS: {class_id}")
            all_results.append((class_id, "SUCCESS"))
        except Exception as e:
            print(f"--> FAILED: {class_id} raised {e!r}")
            all_results.append((class_id, f"FAILED ({e!r})"))
            # Fail the test for pytest reporting
            raise
    else:
        print(f"--> SKIPPED: {class_id} has no marshall_class()")
        all_results.append((class_id, "SKIPPED"))
        pytest.skip(f"{class_id} has no marshall_class()")