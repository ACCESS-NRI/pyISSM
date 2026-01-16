"""
Primary class for all ISSM model interactions.
"""
import numpy as np
import copy
from pyissm.model import classes, mesh, param
from pyissm.tools import wrappers

import sys

class Model():
    """
    ISSM Model Class.

    This class defines a high-level container for all components of an ISSM (Ice Sheet System Model) model.
    It initializes a collection of model components, each of which may store inputs, settings, and results related
    to various aspects of the ice sheet simulation.

    Parameters
    ----------
    None.

    Attributes
    ----------
    mesh : classes.mesh.mesh2d()
        Mesh properties.
    mask : classes.mask.mask2d()
        Defines grounded and floating elements.
    geometry : classes.geometry.geometry2d()
        Surface elevation, bedrock topography, ice thickness, etc.
    constants : classes.constants()
        Physical constants.
    smb : classes.smb.default()
        Surface mass balance.
    basalforcings : classes.basalforcings.default()
        Bed forcings.
    materials : classes.materials.ice()
        Material properties.
    damage : classes.damage()
        Damage propagation laws.
    friction : classes.friction.default()
        Basal friction / drag properties.
    flowequation : classes.flowequation()
        Flow equations.
    timestepping : classes.timestepping.default()
        Timestepping for transient models.
    initialization : classes.initialization()
        Initial guess / state.
    rifts : classes.rifts()
        Rifts properties.
    solidearth : classes.solidearth.earth()
        Solidearth inputs and settings.
    dsl : classes.dsl.default()
        Dynamic sea level.
    debug : classes.debug()
        Debugging tools (valgrind, gprof).
    verbose : classes.verbose()
        Verbosity level in solve.
    settings : classes.issmsettings()
        Settings properties.
    toolkits : None
        PETSc options for each solution.
    cluster : None
        Cluster parameters (number of CPUs, etc.).
    balancethickness : classes.balancethickness()
        Parameters for balancethickness solution.
    stressbalance : classes.stressbalance()
        Parameters for stressbalance solution.
    groundingline : classes.groundingline()
        Parameters for groundingline solution.
    hydrology : classes.hydrology.shreve()
        Parameters for hydrology solution.
    masstransport : classes.masstransport()
        Parameters for masstransport solution.
    thermal : classes.thermal()
        Parameters for thermal solution.
    steadystate : classes.steadystate()
        Parameters for steadystate solution.
    transient : classes.transient()
        Parameters for transient solution.
    levelset : classes.levelset()
        Parameters for moving boundaries (level-set method).
    calving : classes.calving.default()
        Parameters for calving.
    frontalforcings : classes.frontalforcings.default()
        Parameters for frontalforcings.
    love : classes.love.default()
        Parameters for love solution.
    esa : classes.esa()
        Parameters for elastic adjustment solution.
    sampling : classes.sampling()
        Parameters for stochastic sampler.
    autodiff : classes.autodiff()
        Automatic differentiation parameters.
    inversion : classes.inversion.default()
        Parameters for inverse methods.
    qmu : classes.qmu.default()
        Dakota properties.
    amr : classes.amr()
        Adaptive mesh refinement properties.
    results : classes.results.default()
        Model results.
    outputdefinition : classes.outputdefinition()
        Output definition.
    radaroverlay : classes.radaroverlay()
        Radar image for plot overlay.
    miscellaneous : classes.miscellaneous()
        Miscellaneous fields.
    stochasticforcing : classes.stochasticforcing()
        Stochasticity applied to model forcings.
    """

    def __init__(self):

        ## Initialise all as None
        self.mesh = classes.mesh.mesh2d()
        self.mask = classes.mask()
        self.geometry = classes.geometry()
        self.constants = classes.constants()
        self.smb = classes.smb.default()
        self.basalforcings = classes.basalforcings.default()
        self.materials = classes.materials.ice()
        self.damage = classes.damage()
        self.friction = classes.friction.default()
        self.flowequation = classes.flowequation()
        self.timestepping = classes.timestepping.default()
        self.initialization = classes.initialization()
        self.rifts = classes.rifts()
        self.dsl = classes.dsl.default()
        self.solidearth = classes.solidearth.earth()
        self.debug = classes.debug()
        self.verbose = classes.verbose()
        self.settings = classes.issmsettings()
        self.toolkits = classes.toolkits()
        self.cluster = classes.cluster.generic()
        self.balancethickness = classes.balancethickness()
        self.stressbalance = classes.stressbalance()
        self.groundingline = classes.groundingline()
        self.hydrology = classes.hydrology.shreve()
        self.debris = classes.debris()
        self.masstransport = classes.masstransport()
        self.thermal = classes.thermal()
        self.steadystate = classes.steadystate()
        self.transient = classes.transient()
        self.levelset = classes.levelset()
        self.calving = classes.calving.default()
        self.frontalforcings = classes.frontalforcings.default()
        self.love = classes.love.default()
        self.esa = classes.esa()
        self.sampling = classes.sampling()
        self.autodiff = classes.autodiff()
        self.inversion = classes.inversion.default()
        self.qmu = classes.qmu.default()
        self.amr = classes.amr()
        self.results = classes.results.default()
        self.outputdefinition = classes.outputdefinition()
        self.radaroverlay = classes.radaroverlay()
        self.miscellaneous = classes.miscellaneous()
        self.private = classes.private()
        self.stochasticforcing = classes.stochasticforcing()

    # Define repr
    def __repr__(self):
        # Largely consistent with current MATLAB setup
        s = '%19s %-23s %s' % ('ISSM Model Class', '', '')
        s = '%s\n%s' % (s, '%19s %-23s %s' % ('', '', ''))
        s = '%s\n%s' % (s, '%19s:  %-23s %s' % ('mesh', 'mesh properties', ''))
        s = '%s\n%s' % (s, '%19s:  %-23s %s' % ('mask', 'defines grounded and gloating elements', ''))
        s = '%s\n%s' % (s, '%19s:  %-23s %s' % ('geometry', 'surface elevation, bedrock topography, ice thickness, ...', ''))
        s = '%s\n%s' % (s, '%19s:  %-23s %s' % ('constants', 'physical constants', ''))
        s = '%s\n%s' % (s, '%19s:  %-23s %s' % ('smb', 'surface mass balance', ''))
        s = '%s\n%s' % (s, '%19s:  %-23s %s' % ('basalforcings', 'bed forcings', ''))
        s = '%s\n%s' % (s, '%19s:  %-23s %s' % ('materials', 'material properties', ''))
        s = '%s\n%s' % (s, '%19s:  %-23s %s' % ('damage', 'damage propagation laws', ''))
        s = '%s\n%s' % (s, '%19s:  %-23s %s' % ('friction', 'basal friction / drag properties', ''))
        s = '%s\n%s' % (s, '%19s:  %-23s %s' % ('flowequation', 'flow equations', ''))
        s = '%s\n%s' % (s, '%19s:  %-23s %s' % ('timestepping', 'timestepping for transient models', ''))
        s = '%s\n%s' % (s, '%19s:  %-23s %s' % ('initialization', 'initial guess / state', ''))
        s = '%s\n%s' % (s, '%19s:  %-23s %s' % ('rifts', 'rifts properties', ''))
        s = '%s\n%s' % (s, '%19s:  %-23s %s' % ('solidearth', 'solidearth inputs and settings', ''))
        s = '%s\n%s' % (s, '%19s:  %-23s %s' % ('dsl', 'dynamic sea level', ''))
        s = '%s\n%s' % (s, '%19s:  %-23s %s' % ('debug', 'debugging tools (valgrind, gprof', ''))
        s = '%s\n%s' % (s, '%19s:  %-23s %s' % ('verbose', 'verbosity level in solve', ''))
        s = '%s\n%s' % (s, '%19s:  %-23s %s' % ('settings', 'settings properties', ''))
        s = '%s\n%s' % (s, '%19s:  %-23s %s' % ('toolkits', 'PETSc options for each solution', ''))
        s = '%s\n%s' % (s, '%19s:  %-23s %s' % ('cluster', 'cluster parameters (number of CPUs...)', ''))
        s = '%s\n%s' % (s, '%19s:  %-23s %s' % ('balancethickness', 'parameters for balancethickness solution', ''))
        s = '%s\n%s' % (s, '%19s:  %-23s %s' % ('stressbalance', 'parameters for stressbalance solution', ''))
        s = '%s\n%s' % (s, '%19s:  %-23s %s' % ('groundingline', 'parameters for groundingline solution', ''))
        s = '%s\n%s' % (s, '%19s:  %-23s %s' % ('hydrology', 'parameters for hydrology solution', ''))
        s = '%s\n%s' % (s, '%19s:  %-23s %s' % ('masstransport', 'parameters for masstransport solution', ''))
        s = '%s\n%s' % (s, '%19s:  %-23s %s' % ('thermal', 'parameters fo thermal solution', ''))
        s = '%s\n%s' % (s, '%19s:  %-23s %s' % ('steadystate', 'parameters for steadystate solution', ''))
        s = '%s\n%s' % (s, '%19s:  %-23s %s' % ('transient', 'parameters for transient solution', ''))
        s = '%s\n%s' % (s, '%19s:  %-23s %s' % ('levelset', 'parameters for moving boundaries (level-set method)', ''))
        s = '%s\n%s' % (s, '%19s:  %-23s %s' % ('calving', 'parameters for calving', ''))
        s = '%s\n%s' % (s, '%19s:  %-23s %s' % ('frontalforcings', 'parameters for frontalforcings', ''))
        s = '%s\n%s' % (s, '%19s:  %-23s %s' % ('esa', 'parameters for elastic adjustment solution', ''))
        s = '%s\n%s' % (s, '%19s:  %-23s %s' % ('sampling', 'parameters for stochastic sampler', ''))
        s = '%s\n%s' % (s, '%19s:  %-23s %s' % ('love', 'parameters for love solution', ''))
        s = '%s\n%s' % (s, '%19s:  %-23s %s' % ('autodiff', 'automatic differentiation parameters', ''))
        s = '%s\n%s' % (s, '%19s:  %-23s %s' % ('inversion', 'parameters for inverse methods', ''))
        s = '%s\n%s' % (s, '%19s:  %-23s %s' % ('qmu', 'Dakota properties', ''))
        s = '%s\n%s' % (s, '%19s:  %-23s %s' % ('amr', 'adaptive mesh refinement properties', ''))
        s = '%s\n%s' % (s, '%19s:  %-23s %s' % ('outputdefinition', 'output definition', ''))
        s = '%s\n%s' % (s, '%19s:  %-23s %s' % ('results', 'modelresults', ''))
        s = '%s\n%s' % (s, '%19s:  %-23s %s' % ('radaroverlay', 'radar image for plot overlay', ''))
        s = '%s\n%s' % (s, '%19s:  %-23s %s' % ('miscellaneous', 'miscellaneous fields', ''))
        s = '%s\n%s' % (s, '%19s:  %-23s %s' % ('stochasticforcing', 'stochasticity applied to model forcings', ''))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM Model Class'
        return s
    
    def check_message(self, string):
        """
        Notify about a model consistency error, update internal state, and return the instance.

        This method prints a formatted consistency error message to standard output,
        marks the instance as inconsistent by setting ``self.private.isconsistent``
        to ``False``, and returns the instance to allow for method chaining.

        Parameters
        ----------
        string : str
            Human-readable description of the consistency error. This will be inserted
            into the printed message: ``Model consistency error: {string}``.

        Returns
        -------
        self
            The same instance on which the method was called, enabling fluent/chained
            calls.

        Notes
        -----
        This method has the side effect of mutating the instance state (``self.private.isconsistent``),
        and it performs output via ``print``. It does not raise exceptions.

        Examples
        --------
        >>> obj.check_message("missing parameter")
        Model consistency error: missing parameter
        >>> obj.private.isconsistent
        False
        """
        print(f'Model consistency error: {string}')
        self.private.isconsistent = False
        return self
    
    def model_class_names(self):
        """
        Return a sorted list of registered model class attribute names.

        The method inspects the instance attributes and returns those whose
        classes are registered in ``classes.class_registry.CLASS_REGISTRY``.

        Returns
        -------
        list of str
            Sorted list of attribute names corresponding to registered model classes.
        """
        registered_classes = set(classes.class_registry.CLASS_REGISTRY.values())
        names = [
            name for name, obj in vars(self).items()
            if obj.__class__ in registered_classes
        ]
        return sorted(names)

    # Define state
    def __getstate__(self):
        return self.__dict__.copy()
    
    # Extract a portion of existing model
    def extract(self, area):
        """
        Extract a submodel from a larger model based on a domain or flag list.

        This routine extracts a submodel from a bigger model with respect to a given 
        contour. The contour must be followed by the corresponding .exp file or a 
        flags list. It can either be a domain file (argus type, .exp extension), or 
        an array of element flags. If the user wants every element outside the domain 
        to be extracted, add '~' to the name of the domain file (e.g., '~HO.exp'). 
        An empty string '' will be considered as an empty domain. A string 'all' will 
        be considered as the entire domain.

        The function performs the following operations:
        - Flags elements inside the specified area
        - Removes elements with all three nodes in excluded regions
        - Renumbers elements and nodes to maintain consistency
        - Updates mesh connectivity and vertex information
        - Adjusts boundary conditions at extraction boundaries
        - Handles 2D/3D mesh types and their specific properties
        - Processes results and output definitions if present

        Parameters
        ----------
        area : str or array_like
            Domain specification. Can be:
            - A domain file path (argus type, .exp extension)
            - A domain file path prefixed with '~' to invert the domain
            - An array of element flags (boolean or integer)
            - An empty string '' for empty domain
            - The string 'all' for entire domain

        Returns
        -------
        md2 : Model
            Extracted submodel containing only the elements and nodes within the 
            specified area. The model includes updated mesh properties, renumbered 
            elements and nodes, adjusted boundary conditions, and extracted results 
            and output definitions if present.

        Raises
        ------
        RuntimeError
            If the extracted model is empty (no elements found in the specified area).

        See Also
        --------
        extrude : Extrude model in vertical direction
        collapse : Collapse model layers

        Examples
        --------
        >>> md2 = extract(md, 'Domain.exp')
        >>> md3 = extract(md, '~Domain.exp')  # Extract outside domain
        >>> md4 = extract(md, flag_array)     # Extract based on flag array
        """

        ## NOTE: This function is taken directly from $ISSM_DIR/src/m/classes/model.py with only minor modifications for pyISSM integration.

        # Copy model
        md1 = copy.deepcopy(self)

        # Get elements that are inside area
        flag_elem = mesh.flag_elements(md1, area)
        if not np.any(flag_elem):
            raise RuntimeError('extracted model is empty')

        # Kick out all elements with 3 dirichlets
        spc_elem = np.nonzero(np.logical_not(flag_elem))[0]
        spc_node = np.unique(md1.mesh.elements[spc_elem, :]) - 1
        flag = np.ones(md1.mesh.numberofvertices)
        flag[spc_node] = 0
        pos = np.nonzero(np.logical_not(np.sum(flag[md1.mesh.elements - 1], axis=1)))[0]
        flag_elem[pos] = 0

        # Extracted elements and nodes lists
        pos_elem = np.nonzero(flag_elem)[0]
        pos_node = np.unique(md1.mesh.elements[pos_elem, :]) - 1

        # Keep track of some fields
        numberofvertices1 = md1.mesh.numberofvertices
        numberofelements1 = md1.mesh.numberofelements
        numberofvertices2 = np.size(pos_node)
        numberofelements2 = np.size(pos_elem)
        flag_node = np.zeros(numberofvertices1)
        flag_node[pos_node] = 1

        # Create Pelem and Pnode (transform old nodes in new nodes and same thing for the elements)
        Pelem = np.zeros(numberofelements1, int)
        Pelem[pos_elem] = np.arange(1, numberofelements2 + 1)
        Pnode = np.zeros(numberofvertices1, int)
        Pnode[pos_node] = np.arange(1, numberofvertices2 + 1)

        # Renumber the elements (some nodes won't exist anymore)
        elements_1 = copy.deepcopy(md1.mesh.elements)
        elements_2 = elements_1[pos_elem, :]
        elements_2[:, 0] = Pnode[elements_2[:, 0] - 1]
        elements_2[:, 1] = Pnode[elements_2[:, 1] - 1]
        elements_2[:, 2] = Pnode[elements_2[:, 2] - 1]
        if md1.mesh.__class__.__name__ == 'mesh3dprisms':
            elements_2[:, 3] = Pnode[elements_2[:, 3] - 1]
            elements_2[:, 4] = Pnode[elements_2[:, 4] - 1]
            elements_2[:, 5] = Pnode[elements_2[:, 5] - 1]

        # Ok, now create the new model
        # Take every field from model
        md2 = copy.deepcopy(md1)

        # Automatically modify fields
        # Loop over model fields
        md_fieldnames = vars(md1)
        for md_fieldname in md_fieldnames:
            # Get field
            field = getattr(md1, md_fieldname)
            fieldsize = np.shape(field)
            if hasattr(field, '__dict__') and md_fieldname not in ['results']: # recursive call
                obj_fieldnames = vars(field)
                for obj_fieldname in obj_fieldnames:
                    # Get field
                    field = getattr(getattr(md1, md_fieldname), obj_fieldname)
                    fieldsize = np.shape(field)
                    if len(fieldsize):
                        # size = number of nodes * n
                        if fieldsize[0] == numberofvertices1:
                            setattr(getattr(md2, md_fieldname), obj_fieldname, field[pos_node])
                        elif fieldsize[0] == numberofvertices1 + 1:
                            setattr(getattr(md2, md_fieldname), obj_fieldname, np.vstack((field[pos_node], field[-1, :])))
                        # size = number of elements * n
                        elif fieldsize[0] == numberofelements1:
                            setattr(getattr(md2, md_fieldname), obj_fieldname, field[pos_elem])
            else:
                if len(fieldsize):
                    # size = number of nodes * n
                    if fieldsize[0] == numberofvertices1:
                        setattr(md2, md_fieldname, field[pos_node])
                    elif fieldsize[0] == numberofvertices1 + 1:
                        setattr(md2, md_fieldname, np.hstack((field[pos_node], field[-1, :])))
                    # size = number of elements * n
                    elif fieldsize[0] == numberofelements1:
                        setattr(md2, md_fieldname, field[pos_elem])

        # Modify some specific fields
        # mesh
        md2.mesh.numberofelements = numberofelements2
        md2.mesh.numberofvertices = numberofvertices2
        md2.mesh.elements = elements_2

        # mesh.uppervertex mesh.lowervertex
        if isinstance(md1.mesh, classes.mesh.mesh3dprisms):
            md2.mesh.uppervertex = md1.mesh.uppervertex[pos_node]
            pos = np.where(~np.isnan(md2.mesh.uppervertex))[0]
            md2.mesh.uppervertex[pos] = Pnode[md2.mesh.uppervertex[pos].astype(int) - 1]

            md2.mesh.lowervertex = md1.mesh.lowervertex[pos_node]
            pos = np.where(~np.isnan(md2.mesh.lowervertex))[0]
            md2.mesh.lowervertex[pos] = Pnode[md2.mesh.lowervertex[pos].astype(int) - 1]

            md2.mesh.upperelements = md1.mesh.upperelements[pos_elem]
            pos = np.where(~np.isnan(md2.mesh.upperelements))[0]
            md2.mesh.upperelements[pos] = Pelem[md2.mesh.upperelements[pos].astype(int) - 1]

            md2.mesh.lowerelements = md1.mesh.lowerelements[pos_elem]
            pos = np.where(~np.isnan(md2.mesh.lowerelements))[0]
            md2.mesh.lowerelements[pos] = Pelem[md2.mesh.lowerelements[pos].astype(int) - 1]

        # Initial 2d mesh
        if isinstance(md1.mesh, classes.mesh.mesh3dprisms):
            flag_elem_2d = flag_elem[np.arange(0, md1.mesh.numberofelements2d)]
            pos_elem_2d = np.nonzero(flag_elem_2d)[0]
            flag_node_2d = flag_node[np.arange(0, md1.mesh.numberofvertices2d)]
            pos_node_2d = np.nonzero(flag_node_2d)[0]

            md2.mesh.numberofelements2d = np.size(pos_elem_2d)
            md2.mesh.numberofvertices2d = np.size(pos_node_2d)
            md2.mesh.elements2d = md1.mesh.elements2d[pos_elem_2d, :]
            md2.mesh.elements2d[:, 0] = Pnode[md2.mesh.elements2d[:, 0] - 1]
            md2.mesh.elements2d[:, 1] = Pnode[md2.mesh.elements2d[:, 1] - 1]
            md2.mesh.elements2d[:, 2] = Pnode[md2.mesh.elements2d[:, 2] - 1]

            md2.mesh.x2d = md1.mesh.x[pos_node_2d]
            md2.mesh.y2d = md1.mesh.y[pos_node_2d]

        # Edges
        if md1.mesh.domain_type() == '2Dhorizontal':
            if np.ndim(md2.mesh.edges) > 1 and np.size(md2.mesh.edges, axis=1) > 1: # do not use ~isnan because there are some np.nans...
                # Renumber first two columns
                pos = np.nonzero(md2.mesh.edges[:, 3] != -1)[0]
                md2.mesh.edges[:, 0] = Pnode[md2.mesh.edges[:, 0] - 1]
                md2.mesh.edges[:, 1] = Pnode[md2.mesh.edges[:, 1] - 1]
                md2.mesh.edges[:, 2] = Pelem[md2.mesh.edges[:, 2] - 1]
                md2.mesh.edges[pos, 3] = Pelem[md2.mesh.edges[pos, 3] - 1]
                # Remove edges when the 2 vertices are not in the domain
                md2.mesh.edges = md2.mesh.edges[np.nonzero(np.logical_and(md2.mesh.edges[:, 0], md2.mesh.edges[:, 1]))[0], :]
                # Replace all zeros by - 1 in the last two columns
                pos = np.nonzero(md2.mesh.edges[:, 2] == 0)[0]
                md2.mesh.edges[pos, 2] = -1
                pos = np.nonzero(md2.mesh.edges[:, 3] == 0)[0]
                md2.mesh.edges[pos, 3] = -1
                # Invert - 1 on the third column with last column (also invert first two columns!)
                pos = np.nonzero(md2.mesh.edges[:, 2] == -1)[0]
                md2.mesh.edges[pos, 2] = md2.mesh.edges[pos, 3]
                md2.mesh.edges[pos, 3] = -1
                values = md2.mesh.edges[pos, 1]
                md2.mesh.edges[pos, 1] = md2.mesh.edges[pos, 0]
                md2.mesh.edges[pos, 0] = values
                # Finally remove edges that do not belong to any element
                pos = np.nonzero(np.logical_and(md2.mesh.edges[:, 2] == -1, md2.mesh.edges[:, 3] == -1))[0]
                md2.mesh.edges = np.delete(md2.mesh.edges, pos, axis=0)

        # Penalties
        if np.any(np.logical_not(np.isnan(md2.stressbalance.vertex_pairing))):
            for i in range(np.size(md1.stressbalance.vertex_pairing, axis=0)):
                md2.stressbalance.vertex_pairing[i, :] = Pnode[md1.stressbalance.vertex_pairing[i, :]]
            md2.stressbalance.vertex_pairing = md2.stressbalance.vertex_pairing[np.nonzero(md2.stressbalance.vertex_pairing[:, 0])[0], :]
        if np.any(np.logical_not(np.isnan(md2.masstransport.vertex_pairing))):
            for i in range(np.size(md1.masstransport.vertex_pairing, axis=0)):
                md2.masstransport.vertex_pairing[i, :] = Pnode[md1.masstransport.vertex_pairing[i, :]]
            md2.masstransport.vertex_pairing = md2.masstransport.vertex_pairing[np.nonzero(md2.masstransport.vertex_pairing[:, 0])[0], :]

        # Recreate segments
        if isinstance(md1.mesh, classes.mesh.mesh2d):
            md2.mesh.vertexconnectivity = wrappers.NodeConnectivity(md2.mesh.elements, md2.mesh.numberofvertices)
            md2.mesh.elementconnectivity = wrappers.ElementConnectivity(md2.mesh.elements, md2.mesh.vertexconnectivity)
            md2.mesh.segments = param.contour_envelope(md2.mesh)
            md2.mesh.vertexonboundary = np.zeros(numberofvertices2, int)
            md2.mesh.vertexonboundary[md2.mesh.segments[:, 0:2] - 1] = 1
        else:
            # First do the connectivity for the contourenvelope in 2d
            md2.mesh.vertexconnectivity = wrappers.NodeConnectivity(md2.mesh.elements2d, md2.mesh.numberofvertices2d)
            md2.mesh.elementconnectivity = wrappers.ElementConnectivity(md2.mesh.elements2d, md2.mesh.vertexconnectivity)
            segments = param.contour_envelope(md2.mesh)
            md2.mesh.vertexonboundary = np.zeros(int(numberofvertices2 / md2.mesh.numberoflayers), int)
            md2.mesh.vertexonboundary[segments[:, 0:2] - 1] = 1
            md2.mesh.vertexonboundary = np.tile(md2.mesh.vertexonboundary, md2.mesh.numberoflayers)
            # Then do it for 3d as usual
            md2.mesh.vertexconnectivity = wrappers.NodeConnectivity(md2.mesh.elements, md2.mesh.numberofvertices)
            md2.mesh.elementconnectivity = wrappers.ElementConnectivity(md2.mesh.elements, md2.mesh.vertexconnectivity)

        # Boundary conditions: Dirichlets on new boundary
        # Catch the elements that have not been extracted
        orphans_elem = np.nonzero(np.logical_not(flag_elem))[0]
        orphans_node = np.unique(md1.mesh.elements[orphans_elem, :]) - 1
        # Figure out which node are on the boundary between md2 and md1
        nodestoflag1 = np.intersect1d(orphans_node, pos_node)
        nodestoflag2 = Pnode[nodestoflag1].astype(int) - 1
        if np.size(md1.stressbalance.spcvx) > 1 and np.size(md1.stressbalance.spcvy) > 2 and np.size(md1.stressbalance.spcvz) > 2:
            if np.size(md1.inversion.vx_obs) > 1 and np.size(md1.inversion.vy_obs) > 1:
                md2.stressbalance.spcvx[nodestoflag2] = md2.inversion.vx_obs[nodestoflag2]
                md2.stressbalance.spcvy[nodestoflag2] = md2.inversion.vy_obs[nodestoflag2]
            else:
                md2.stressbalance.spcvx[nodestoflag2] = np.nan
                md2.stressbalance.spcvy[nodestoflag2] = np.nan
                print('\n!! extract warning: spc values should be checked !!\n\n')
            # Put 0 for vz
            md2.stressbalance.spcvz[nodestoflag2] = 0
        if np.any(np.logical_not(np.isnan(md1.thermal.spctemperature))):
            md2.thermal.spctemperature[nodestoflag2] = 1

        # Results fields
        if md1.results:
            md2.results = classes.results.default()
            for solutionfield, field in list(md1.results.__dict__.items()):
                if isinstance(field, list):
                    setattr(md2.results, solutionfield, [])
                    # Get time step
                    for i, fieldi in enumerate(field):
                        if isinstance(fieldi, classes.results.default) and fieldi:
                            getattr(md2.results, solutionfield).append(classes.results.default())
                            fieldr = getattr(md2.results, solutionfield)[i]
                            # Get subfields
                            for solutionsubfield, subfield in list(fieldi.__dict__.items()):
                                if np.size(subfield) == numberofvertices1:
                                    setattr(fieldr, solutionsubfield, subfield[pos_node])
                                elif np.size(subfield) == numberofelements1:
                                    setattr(fieldr, solutionsubfield, subfield[pos_elem])
                                else:
                                    setattr(fieldr, solutionsubfield, subfield)
                        else:
                            getattr(md2.results, solutionfield).append(None)
                elif isinstance(field, classes.results.default):
                    setattr(md2.results, solutionfield, classes.results.default())
                    if isinstance(field, classes.results.default) and field:
                        fieldr = getattr(md2.results, solutionfield)
                        # Get subfields
                        for solutionsubfield, subfield in list(field.__dict__.items()):
                            if np.size(subfield) == numberofvertices1:
                                setattr(fieldr, solutionsubfield, subfield[pos_node])
                            elif np.size(subfield) == numberofelements1:
                                setattr(fieldr, solutionsubfield, subfield[pos_elem])
                            else:
                                setattr(fieldr, solutionsubfield, subfield)

        # outputdefinitions fields
        if md1.outputdefinition.definitions:
            for solutionfield, field in list(md1.outputdefinition.__dict__.items()):
                if isinstance(field, list):
                    # Get each definition
                    for i, fieldi in enumerate(field):
                        if fieldi:
                            fieldr = getattr(md2.outputdefinition, solutionfield)[i]
                            # Get subfields
                            for solutionsubfield, subfield in list(fieldi.__dict__.items()):
                                if np.size(subfield) == numberofvertices1:
                                    setattr(fieldr, solutionsubfield, subfield[pos_node])
                                elif np.size(subfield) == numberofelements1:
                                    setattr(fieldr, solutionsubfield, subfield[pos_elem])
                                else:
                                    setattr(fieldr, solutionsubfield, subfield)

        # Keep track of pos_node and pos_elem
        md2.mesh.extractedvertices = pos_node + 1
        md2.mesh.extractedelements = pos_elem + 1

        return md2
    
    # Extrude a 2D mesh to 3D mesh
    def extrude(self,
                num_layers = None,
                extrusion_exponent = None,
                lower_exponent = None,
                upper_exponent = None,
                coefficients = None):
        """
        Vertically extrude a 2D mesh to create a 3D mesh.

        Vertically extrude a 2D mesh and create corresponding 3D prism mesh.
        The vertical distribution can follow a polynomial law, two polynomial laws 
        (one for lower part and one for upper part), or be described by a list of 
        coefficients between 0 and 1.

        Parameters
        ----------
        num_layers : int, optional
            Number of vertical layers to create. Required when using polynomial 
            extrusion (single or dual exponent). Must be at least 2.
        extrusion_exponent : float, optional
            Single polynomial exponent for uniform vertical distribution. Must be >= 0.
            When specified, creates extrusion list using this exponent.
        lower_exponent : float, optional
            Polynomial exponent for the lower part of the mesh. Must be >= 0.
            Used in conjunction with `upper_exponent` for non-uniform vertical distribution.
        upper_exponent : float, optional
            Polynomial exponent for the upper part of the mesh. Must be >= 0.
            Used in conjunction with `lower_exponent` for non-uniform vertical distribution.
        coefficients : array_like, optional
            List of coefficients between 0 and 1 defining custom vertical distribution.
            Automatically includes 0 and 1 if not present. Alternative to polynomial extrusion.

        Returns
        -------
        md : Model
            The same model instance with updated 3D mesh and extruded model fields.
            Original 2D mesh is preserved in ``mesh.x2d``, ``mesh.y2d``, 
            ``mesh.elements2d``, etc.

        Raises
        ------
        TypeError
            If extrusion_exponent or lower_exponent/upper_exponent <= 0.
        TypeError
            If coefficients contain values outside [0, 1] range.
        TypeError
            If num_layers < 2.
        TypeError
            If mesh is already 3D (extrude called more than once).

        Notes
        -----
        The function performs the following operations:
        - Creates 3D prism elements from 2D elements
        - Distributes nodes vertically between bed and surface according to extrusion parameters
        - Maintains vertex and element numbering relationships via ``lowervertex``, 
          ``uppervertex``, ``lowerelements``, and ``upperelements``
        - Projects all 2D model fields to 3D using ``extrude`` methods
        - Updates mesh connectivity and boundary information
        - Preserves the original 2D mesh information for reference

        Examples
        --------
        Single exponent (uniform-like distribution):

        >>> md = md.extrude(num_layers=15, extrusion_exponent=1.3)

        Dual exponents (non-uniform distribution):

        >>> md = md.extrude(num_layers=15, lower_exponent=1.3, upper_exponent=1.2)

        Custom coefficients (specific layer distribution):

        >>> md = md.extrude(coefficients=[0, 0.2, 0.5, 0.7, 0.9, 0.95, 1])
        """

        ## NOTE: This function is taken directly from $ISSM_DIR/src/m/classes/model.py with only minor modifications for pyISSM integration.

        md = copy.deepcopy(self)

        # Error checks and Parse inputs
        if coefficients is not None: # list of coefficients
            clist = coefficients
            if any(clist < 0) or any(clist > 1):
                raise TypeError('pyissm.model.Model.extrude: All coefficients must be between 0 and 1')
            clist.extend([0., 1.])
            clist.sort()
            extrusionlist = list(set(clist))
            numlayers = len(extrusionlist)

        elif extrusion_exponent is not None: # one polynomial law
            if extrusion_exponent <= 0:
                raise TypeError('pyissm.model.Model.extrude: extrusion_exponent must be >= 0')
            numlayers = num_layers
            extrusionlist = (np.arange(0., float(numlayers - 1) + 1., 1.) / float(numlayers - 1))**extrusion_exponent

        elif lower_exponent is not None or upper_exponent is not None: # two polynomial laws
            numlayers = num_layers
            lowerexp = lower_exponent
            upperexp = upper_exponent

            if lower_exponent <= 0 or upper_exponent <= 0:
                raise TypeError('pyissm.model.Model.extrude:: lower_exponent and upper_exponent must be >= 0')

            lowerextrusionlist = (np.arange(0., 1. + 2. / float(numlayers - 1), 2. / float(numlayers - 1)))**lowerexp / 2.
            upperextrusionlist = (np.arange(0., 1. + 2. / float(numlayers - 1), 2. / float(numlayers - 1)))**upperexp / 2.
            extrusionlist = np.unique(np.concatenate((lowerextrusionlist, 1. - upperextrusionlist)))

        if numlayers < 2:
            raise TypeError('pyissm.model.Model.extrude: num_layers should be at least 2')
        
        if isinstance(md.mesh, classes.mesh3dprisms):
            raise TypeError('pyissm.model.Model.extrude: Model mesh is already 3D. Cannot extrude a 3D mesh (extrude cannot be called more than once)')

        # Initialize with 2d mesh
        mesh2d = md.mesh
        md.mesh = classes.mesh3dprisms()

        md.mesh.x = mesh2d.x
        md.mesh.y = mesh2d.y
        md.mesh.elements = mesh2d.elements
        md.mesh.numberofelements = mesh2d.numberofelements
        md.mesh.numberofvertices = mesh2d.numberofvertices

        md.mesh.lat = mesh2d.lat
        md.mesh.long = mesh2d.long
        md.mesh.epsg = mesh2d.epsg
        md.mesh.scale_factor = mesh2d.scale_factor

        md.mesh.vertexonboundary = mesh2d.vertexonboundary
        md.mesh.vertexconnectivity = mesh2d.vertexconnectivity
        md.mesh.elementconnectivity = mesh2d.elementconnectivity
        md.mesh.average_vertex_connectivity = mesh2d.average_vertex_connectivity

        md.mesh.extractedvertices = mesh2d.extractedvertices
        md.mesh.extractedelements = mesh2d.extractedelements

        md.mesh.segments2d = mesh2d.segments

        x3d = np.empty((0))
        y3d = np.empty((0))
        z3d = np.empty((0)) # the lower node is on the bed
        thickness3d = md.geometry.thickness # thickness and bed for these nodes
        bed3d = md.geometry.base

        # Create the new layers
        for i in range(numlayers):
            x3d = np.concatenate((x3d, md.mesh.x))
            y3d = np.concatenate((y3d, md.mesh.y))
            # Nodes are distributed between bed and surface accordingly to the given exponent
            z3d = np.concatenate((z3d, (bed3d + thickness3d * extrusionlist[i]).reshape(-1)))
        number_nodes3d = np.size(x3d) # Number of 3d nodes for the non-extruded part of the mesh

        # Extrude elements
        elements3d = np.empty((0, 6), int)
        for i in range(numlayers - 1):
            elements3d = np.vstack((elements3d, np.hstack((md.mesh.elements + i * md.mesh.numberofvertices,
                                                           md.mesh.elements + (i + 1) * md.mesh.numberofvertices)))) # create the elements of the 3d mesh for the non-extruded part
        number_el3d = np.size(elements3d, axis=0) # number of 3d nodes for the non-extruded part of the mesh

        # Keep a trace of lower and upper nodes
        lowervertex = np.nan * np.ones(number_nodes3d, int)
        uppervertex = np.nan * np.ones(number_nodes3d, int)
        lowervertex[md.mesh.numberofvertices:] = np.arange(1, (numlayers - 1) * md.mesh.numberofvertices + 1)
        uppervertex[:(numlayers - 1) * md.mesh.numberofvertices] = np.arange(md.mesh.numberofvertices + 1, number_nodes3d + 1)
        md.mesh.lowervertex = lowervertex
        md.mesh.uppervertex = uppervertex

        # Same for lower and upper elements
        lowerelements = np.nan * np.ones(number_el3d, int)
        upperelements = np.nan * np.ones(number_el3d, int)
        lowerelements[md.mesh.numberofelements:] = np.arange(1, (numlayers - 2) * md.mesh.numberofelements + 1)
        upperelements[:(numlayers - 2) * md.mesh.numberofelements] = np.arange(md.mesh.numberofelements + 1, (numlayers - 1) * md.mesh.numberofelements + 1)
        md.mesh.lowerelements = lowerelements
        md.mesh.upperelements = upperelements

        # Save old mesh
        md.mesh.x2d = md.mesh.x
        md.mesh.y2d = md.mesh.y
        md.mesh.elements2d = md.mesh.elements
        md.mesh.numberofelements2d = md.mesh.numberofelements
        md.mesh.numberofvertices2d = md.mesh.numberofvertices

        # Build global 3d mesh
        md.mesh.elements = elements3d
        md.mesh.x = x3d
        md.mesh.y = y3d
        md.mesh.z = z3d
        md.mesh.numberofelements = number_el3d
        md.mesh.numberofvertices = number_nodes3d
        md.mesh.numberoflayers = numlayers

        # Ok, now deal with the other fields from the 2d mesh
        # Bed info and surface info
        md.mesh.vertexonbase = mesh.project_3d(md, vector = np.ones(md.mesh.numberofvertices2d, int), type = 'node', layer = 1)
        md.mesh.vertexonsurface = mesh.project_3d(md, vector = np.ones(md.mesh.numberofvertices2d, int), type = 'node', layer = md.mesh.numberoflayers)
        md.mesh.vertexonboundary = mesh.project_3d(md, vector = md.mesh.vertexonboundary, type = 'node')

        # lat/long
        md.mesh.lat = mesh.project_3d(md, vector = md.mesh.lat, type = 'node')
        md.mesh.long = mesh.project_3d(md, vector = md.mesh.long, type = 'node')
        md.mesh.scale_factor = mesh.project_3d(md, vector = md.mesh.scale_factor, type = 'node')

        # Project model fields
        md.geometry.extrude(md)
        md.friction.extrude(md)
        md.inversion.extrude(md)
        md.smb.extrude(md)
        md.initialization.extrude(md)
        md.flowequation.extrude(md)
        md.stressbalance.extrude(md)
        md.thermal.extrude(md)
        md.masstransport.extrude(md)
        md.levelset.extrude(md)
        md.calving.extrude(md)
        md.frontalforcings.extrude(md)
        md.hydrology.extrude(md)
        md.debris.extrude(md)
        md.solidearth.extrude(md)
        md.dsl.extrude(md)
        md.stochasticforcing.extrude(md)
        md.materials.extrude(md)
        md.damage.extrude(md)
        md.mask.extrude(md)
        md.qmu.extrude(md)
        md.basalforcings.extrude(md)
        md.outputdefinition.extrude(md)

        # Update connectivity
        if not np.isnan(md.mesh.elementconnectivity).all():
            ne2d = md.mesh.numberofelements2d

            ## Use floats to allow NaNs
            elemconn = md.mesh.elementconnectivity.astype(float)

            ## Replicate
            elemconn = np.tile(elemconn, (num_layers -1, 1))

            ## Repalce 0 with nan
            elemconn[elemconn == 0] = np.nan

            ## Shift layer numbers:
            for i in range(2, num_layers):
                start = (i - 1) * ne2d
                end = i * ne2d
                elemconn[start:end, :] += ne2d

            ## Replace nan with 0
            elemconn [np.isnan(elemconn)] = 0

            ## Convert back to int
            md.mesh.elementconnectivity = elemconn.astype(md.mesh.elementconnectivity.dtype)

        if md.mesh.average_vertex_connectivity <= 25:
            md.mesh.average_vertex_connectivity = 100

        return md
    
    # Collapse 3D mesh to 2D mesh
    def collapse(self):
        """
        Collapse a 3D mesh into a 2D mesh.

        This routine collapses a 3D model into a 2D model and collapses all the
        fields of the 3D model by taking their depth-averaged values.

        The function performs the following operations:
        - Projects 3D friction coefficients to 2D (at bedrock level)
        - Averages or projects 3D observations and initialization fields to 2D
        - Collapses boundary conditions from 3D to 2D
        - Depth-averages velocity and material properties
        - Projects 3D results and output definitions to 2D
        - Rebuilds mesh connectivity for the 2D mesh

        Returns
        -------
        md : Model
            The same model instance with updated 2D mesh and collapsed fields.
            Original 3D mesh information is preserved.

        Raises
        ------
        Exception
            If the model does not contain a 3D mesh.
        Exception
            If an unsupported friction type is encountered.

        Examples
        --------
        >>> md = md.collapse()
        """

        ## NOTE: This function is taken directly from $ISSM_DIR/src/m/classes/model.py with only minor modifications for pyISSM integration.

        md = copy.deepcopy(self)

        # Check that the model is really a 3d model
        if md.mesh.domain_type().lower() != '3d':
            raise Exception('pyissm.model.Model.collapse: md must contain a 3D mesh')

        # Start with changing all the fields from the 3d mesh

        # Dealing with the friction law
        # Drag is limited to nodes that are on the bedrock.
        if isinstance(md.friction, classes.friction.default):
            md.friction.coefficient = mesh.project_2d(md, md.friction.coefficient, 1)
            md.friction.p = mesh.project_2d(md, md.friction.p, 1)
            md.friction.q = mesh.project_2d(md, md.friction.q, 1)
        elif isinstance(md.friction, classes.friction.coulomb):
            md.friction.coefficient = mesh.project_2d(md, md.friction.coefficient, 1)
            md.friction.coefficientcoulomb = mesh.project_2d(md, md.friction.coefficientcoulomb, 1)
            md.friction.p = mesh.project_2d(md, md.friction.p, 1)
            md.friction.q = mesh.project_2d(md, md.friction.q, 1)
        elif isinstance(md.friction, classes.friction.hydro):
            md.friction.q = mesh.project_2d(md, md.friction.q, 1)
            md.friction.C = mesh.project_2d(md, md.friction.C, 1)
            md.friction.As = mesh.project_2d(md, md.friction.As, 1)
            md.friction.effective_pressure = mesh.project_2d(md, md.friction.effective_pressure, 1)
        elif isinstance(md.friction, classes.friction.waterlayer):
            md.friction.coefficient = mesh.project_2d(md, md.friction.coefficient, 1)
            md.friction.p = mesh.project_2d(md, md.friction.p, 1)
            md.friction.q = mesh.project_2d(md, md.friction.q, 1)
            md.friction.water_layer = mesh.project_2d(md, md.friction.water_layer, 1)
        elif isinstance(md.friction, classes.friction.weertman):
            md.friction.C = mesh.project_2d(md, md.friction.C, 1)
            md.friction.m = mesh.project_2d(md, md.friction.m, 1)
        elif isinstance(md.friction, classes.friction.weertmantemp):
            md.friction.C = mesh.project_2d(md, md.friction.C, 1)
            md.friction.m = mesh.project_2d(md, md.friction.m, 1)
        else:
            raise Exception('pyissm.model.Model.collapse: Friction type not supported for collapse.')

        # Observations
        if not np.isnan(md.inversion.vx_obs).all():
            md.inversion.vx_obs = mesh.project_2d(md, md.inversion.vx_obs, md.mesh.numberoflayers)
        if not np.isnan(md.inversion.vy_obs).all():
            md.inversion.vy_obs = mesh.project_2d(md, md.inversion.vy_obs, md.mesh.numberoflayers)
        if not np.isnan(md.inversion.vel_obs).all():
            md.inversion.vel_obs = mesh.project_2d(md, md.inversion.vel_obs, md.mesh.numberoflayers)
        if not np.isnan(md.inversion.thickness_obs).all():
            md.inversion.thickness_obs = mesh.project_2d(md, md.inversion.thickness_obs, md.mesh.numberoflayers)
        if not np.isnan(md.inversion.cost_functions_coefficients).all():
            md.inversion.cost_functions_coefficients = mesh.project_2d(md, md.inversion.cost_functions_coefficients, md.mesh.numberoflayers)
        if isinstance(md.inversion.min_parameters, np.ndarray) and md.inversion.min_parameters.size > 1:
            md.inversion.min_parameters = mesh.project_2d(md, md.inversion.min_parameters, md.mesh.numberoflayers)
        if isinstance(md.inversion.max_parameters, np.ndarray) and md.inversion.max_parameters.size > 1:
            md.inversion.max_parameters = mesh.project_2d(md, md.inversion.max_parameters, md.mesh.numberoflayers)
        if isinstance(md.smb, classes.smb.default) and not np.isnan(md.smb.mass_balance).all():
            md.smb.mass_balance = mesh.project_2d(md, md.smb.mass_balance, md.mesh.numberoflayers)
        elif isinstance(md.smb, classes.smb.henning) and not np.isnan(md.smb.smbref).all():
            md.smb.smbref = mesh.project_2d(md, md.smb.smbref, md.mesh.numberoflayers)

        # Results
        if not np.isnan(md.initialization.vx).all():
            md.initialization.vx = mesh.depth_average(md, md.initialization.vx)
        if not np.isnan(md.initialization.vy).all():
            md.initialization.vy = mesh.depth_average(md, md.initialization.vy)
        if not np.isnan(md.initialization.vz).all():
            md.initialization.vz = mesh.depth_average(md, md.initialization.vz)
        if not np.isnan(md.initialization.vel).all():
            md.initialization.vel = mesh.depth_average(md, md.initialization.vel)
        if not np.isnan(md.initialization.temperature).all():
            md.initialization.temperature = mesh.depth_average(md, md.initialization.temperature)
        if not np.isnan(md.initialization.pressure).all():
            md.initialization.pressure = mesh.project_2d(md, md.initialization.pressure, 1)
        if not np.isnan(md.initialization.sediment_head).all():
            md.initialization.sediment_head = mesh.project_2d(md, md.initialization.sediment_head, 1)
        if not np.isnan(md.initialization.epl_head).all():
            md.initialization.epl_head = mesh.project_2d(md, md.initialization.epl_head, 1)
        if not np.isnan(md.initialization.epl_thickness).all():
            md.initialization.epl_thickness = mesh.project_2d(md, md.initialization.epl_thickness, 1)
        if not np.isnan(md.initialization.waterfraction).all():
            md.initialization.waterfraction = mesh.project_2d(md, md.initialization.waterfraction, 1)
        if not np.isnan(md.initialization.watercolumn).all():
            md.initialization.watercolumn = mesh.project_2d(md, md.initialization.watercolumn, 1)
        if not np.isnan(md.initialization.debris).all():
            md.initialization.debris = mesh.project_2d(md, md.initialization.debris, 1)

        # Element types
        if not np.isnan(md.flowequation.element_equation).all():
            md.flowequation.element_equation = mesh.project_2d(md, md.flowequation.element_equation, 1)
            md.flowequation.vertex_equation = mesh.project_2d(md, md.flowequation.vertex_equation, 1)
            md.flowequation.borderSSA = mesh.project_2d(md, md.flowequation.borderSSA, 1)
            md.flowequation.borderHO = mesh.project_2d(md, md.flowequation.borderHO, 1)
            md.flowequation.borderFS = mesh.project_2d(md, md.flowequation.borderFS, 1)

        # Boundary conditions
        md.stressbalance.spcvx = mesh.project_2d(md, md.stressbalance.spcvx, md.mesh.numberoflayers)
        md.stressbalance.spcvy = mesh.project_2d(md, md.stressbalance.spcvy, md.mesh.numberoflayers)
        md.stressbalance.spcvz = mesh.project_2d(md, md.stressbalance.spcvz, md.mesh.numberoflayers)
        md.stressbalance.referential = mesh.project_2d(md, md.stressbalance.referential, md.mesh.numberoflayers)
        md.stressbalance.loadingforce = mesh.project_2d(md, md.stressbalance.loadingforce, md.mesh.numberoflayers)

        # TODO:
        # - Check if md.mesh.numberoflayershould really be offset by 1.
        # - Find out why md.masstransport.spcthickness is not offset, but the
        #   other fields are.
        # - If offset is required, figure out if it can be abstracted away to
        #   another part of the API.
        #
        if np.size(md.masstransport.spcthickness) > 1:
            md.masstransport.spcthickness = mesh.project_2d(md, md.masstransport.spcthickness, md.mesh.numberoflayers)
        if np.size(md.damage.spcdamage) > 1:  # and not np.isnan(md.damage.spcdamage).all():
            md.damage.spcdamage = mesh.project_2d(md, md.damage.spcdamage, md.mesh.numberoflayers - 1)
        if np.size(md.levelset.spclevelset) > 1:
            md.levelset.spclevelset = mesh.project_2d(md, md.levelset.spclevelset, md.mesh.numberoflayers - 1)
        md.thermal.spctemperature = mesh.project_2d(md, md.thermal.spctemperature, md.mesh.numberoflayers - 1)

        # hydrologydc variables
        if isinstance(md.hydrology, classes.hydrology.dc):
            # md.hydrology.spcsediment_head = project2d(md, md.hydrology.spcsediment_head, 1)
            # md.hydrology.mask_eplactive_node = project2d(md, md.hydrology.mask_eplactive_node, 1)
            # md.hydrology.sediment_transmitivity = project2d(md, md.hydrology.sediment_transmitivity, 1)
            # md.hydrology.basal_moulin_input = project2d(md, md.hydrology.basal_moulin_input, 1)
            # if md.hydrology.isefficientlayer == 1:
            #     md.hydrology.spcepl_head = project2d(md, md.hydrology.spcepl_head, 1)
            hydrofields = md.hydrology.__dict__.keys()
            for field in hydrofields:
                try:
                    isvector = np.logical_or(np.shape(md.hydrology.__dict__[field])[0] == md.mesh.numberofelements,
                                             np.shape(md.hydrology.__dict__[field])[0] == md.mesh.numberofvertices)
                except IndexError:
                    isvector = False
                #we collapse only fields that are vertices or element based
                if isvector:
                    md.hydrology.__dict__[field] = mesh.project_2d(md, md.hydrology.__dict__[field], 1)

        # materials
        md.materials.rheology_B = mesh.depth_average(md, md.materials.rheology_B)
        md.materials.rheology_n = mesh.project_2d(md, md.materials.rheology_n, 1)

        # dsl
        if np.size(md.dsl.sea_surface_height_above_geoid) > 1:
            md.dsl.sea_surface_height_above_geoid = mesh.project_2d(md, md.dsl.sea_surface_height_above_geoid, 1)
        if np.size(md.dsl.sea_water_pressure_at_sea_floor) > 1:
            md.dsl.sea_water_pressure_at_sea_floor = mesh.project_2d(md, md.dsl.sea_water_pressure_at_sea_floor, 1)

        # damage
        if md.damage.isdamage:
            md.damage.D = mesh.depth_average(md, md.damage.D)

        # Special for thermal modeling
        if not np.isnan(md.basalforcings.groundedice_melting_rate).all():
            md.basalforcings.groundedice_melting_rate = mesh.project_2d(md, md.basalforcings.groundedice_melting_rate, 1)
        if hasattr(md.basalforcings, 'floatingice_melting_rate') and not np.isnan(md.basalforcings.floatingice_melting_rate).all():
            md.basalforcings.floatingice_melting_rate = mesh.project_2d(md, md.basalforcings.floatingice_melting_rate, 1)
        md.basalforcings.geothermalflux = mesh.project_2d(md, md.basalforcings.geothermalflux, 1) # bedrock only gets geothermal flux

        if hasattr(md.calving, 'coeff') and not np.isnan(md.calving.coeff).all():
            md.calving.coeff = mesh.project_2d(md, md.calving.coeff, 1)
        if hasattr(md.frontalforcings, 'meltingrate') and not np.isnan(md.frontalforcings.meltingrate).all():
            md.frontalforcings.meltingrate = mesh.project_2d(md, md.frontalforcings.meltingrate, 1)

        # Update of connectivity matrix
        md.mesh.average_vertex_connectivity = 25

        # Collapse the mesh
        nodes2d = md.mesh.numberofvertices2d
        elements2d = md.mesh.numberofelements2d

        # Parameters
        md.geometry.surface = mesh.project_2d(md, md.geometry.surface, 1)
        md.geometry.thickness = mesh.project_2d(md, md.geometry.thickness, 1)
        md.geometry.base = mesh.project_2d(md, md.geometry.base, 1)
        if not np.isnan(md.geometry.bed).all():
            md.geometry.bed = mesh.project_2d(md, md.geometry.bed, 1)
        if not np.isnan(md.mask.ocean_levelset).all():
            md.mask.ocean_levelset = mesh.project_2d(md, md.mask.ocean_levelset, 1)
        if not np.isnan(md.mask.ice_levelset).all():
            md.mask.ice_levelset = mesh.project_2d(md, md.mask.ice_levelset, 1)

        # lat/long
        if np.size(md.mesh.lat) == md.mesh.numberofvertices:
            md.mesh.lat = mesh.project_2d(md, md.mesh.lat, 1)
        if np.size(md.mesh.long) == md.mesh.numberofvertices:
            md.mesh.long = mesh.project_2d(md, md.mesh.long, 1)

        # outputdefinitions
        if md.outputdefinition.definitions:
            for solutionfield, field in list(md.outputdefinition.__dict__.items()):
                if isinstance(field, list):
                    # Get each definition
                    for i, fieldi in enumerate(field):
                        if fieldi:
                            fieldr = getattr(md.outputdefinition, solutionfield)[i]
                            # Get subfields
                            for solutionsubfield, subfield in list(fieldi.__dict__.items()):
                                if np.size(subfield) == md.mesh.numberofvertices:
                                    setattr(fieldr, solutionsubfield, mesh.project_2d(md, subfield, 1))
                                elif np.size(subfield) == md.mesh.numberofelements:
                                    setattr(fieldr, solutionsubfield, mesh.project_2d(md, subfield, 1))

        # Initialize 2d mesh
        mesh2d = classes.mesh.mesh2d()
        mesh2d.x = md.mesh.x2d
        mesh2d.y = md.mesh.y2d
        mesh2d.numberofvertices = md.mesh.numberofvertices2d
        mesh2d.numberofelements = md.mesh.numberofelements2d
        mesh2d.elements = md.mesh.elements2d
        # if not np.isnan(md.mesh.vertexonboundary).all():
        #     mesh.vertexonboundary = project2d(md, md.mesh.vertexonboundary, 1)
        # if not np.isnan(md.mesh.elementconnectivity).all():
        #     mesh.elementconnectivity = project2d(md, md.mesh.elementconnectivity, 1)
        if np.size(md.mesh.lat) == md.mesh.numberofvertices:
            mesh2d.lat = mesh.project_2d(md, md.mesh.lat, 1)
        if np.size(md.mesh.long) == md.mesh.numberofvertices:
            mesh2d.long = mesh.project_2d(md, md.mesh.long, 1)
        mesh.epsg = md.mesh.epsg
        if np.size(md.mesh.scale_factor) == md.mesh.numberofvertices:
            mesh2d.scale_factor = mesh.project_2d(md, md.mesh.scale_factor, 1)
        if hasattr(md.mesh, 'vertexonboundary') and not np.isnan(md.mesh.vertexonboundary).all():
            mesh2d.vertexonboundary = mesh.project_2d(md, md.mesh.vertexonboundary, 1)
        if hasattr(md.mesh, 'elementonboundary') and not np.isnan(md.mesh.elementonboundary).all():
            mesh2d.elementonboundary = mesh.project_2d(md, md.mesh.elementonboundary, 1)
        md.mesh = mesh2d
        md.mesh.vertexconnectivity = wrappers.NodeConnectivity(md.mesh.elements, md.mesh.numberofvertices)
        md.mesh.elementconnectivity = wrappers.ElementConnectivity(md.mesh.elements, md.mesh.vertexconnectivity)
        md.mesh.segments = param.contour_envelope(md.mesh)

        return md