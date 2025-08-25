import numpy as np
from . import param_utils
from . import class_registry
from .. import execute

@class_registry.register_class
class flowequation(class_registry.manage_state):
    """
    Flow equation parameters class for ISSM.

    This class encapsulates parameters for configuring flow equations in the ISSM (Ice Sheet System Model) framework.
    It allows users to select and configure different ice flow approximations such as SIA, SSA, Higher-Order, and Full-Stokes,
    along with their associated numerical methods and boundary conditions.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    isSIA : int, default=0
        Is the Shallow Ice Approximation (SIA) used?
    isL1L2 : int, default=0
        Are L1L2 equations used?
    isSSA : int, default=0
        Is the Shelfy-Stream Approximation (SSA) used?
    isMOLHO : int, default=0
        Are MOno-layer Higher-Order (MOLHO) equations used?
    isHO : int, default=0
        Is the Higher-Order (HO) approximation used?
    isFS : int, default=0
        Are the Full-Stokes (FS) equations used?
    isNitscheBC : int, default=0
        Is weakly imposed condition used?
    FSNitscheGamma : float, default=1e6
        Gamma value for the Nitsche term.
    fe_SSA : str, default='P1'
        Finite Element for SSA: 'P1', 'P1bubble', 'P1bubblecondensed', 'P2'.
    fe_HO : str, default='P1'
        Finite Element for HO: 'P1', 'P1bubble', 'P1bubblecondensed', 'P2'.
    fe_FS : str, default='MINIcondensed'
        Finite Element for FS: 'MINI', 'MINIcondensed', 'TaylorHood', 'XTaylorHood', 'OneLayerP4z', 'CrouzeixRaviart'.
    augmented_lagrangian_r : float, default=1
        Augmented Lagrangian parameter r.
    augmented_lagrangian_rhop : float, default=1
        Augmented Lagrangian parameter rhop.
    augmented_lagrangian_rlambda : float, default=1
        Augmented Lagrangian parameter rlambda.
    augmented_lagrangian_rholambda : float, default=1
        Augmented Lagrangian parameter rholambda.
    XTH_theta : float, default=0
        XTH theta parameter.
    vertex_equation : float, default=nan
        Vertex equation parameter.
    element_equation : float, default=nan
        Element equation parameter.
    borderSSA : float, default=nan
        Border parameter for SSA.
    borderHO : float, default=nan
        Border parameter for HO.
    borderFS : float, default=nan
        Border parameter for FS.

    Methods
    -------
    __init__(self, other=None)
        Initializes the flowequation parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the flowequation parameters.
    __str__(self)
        Returns a short string identifying the class.
    marshall_class(self, fid, prefix, md=None)
        Marshall parameters to a binary file.

    Examples
    --------
    md.flowequation = pyissm.param.flowequation()
    md.flowequation.isSSA = 1
    md.flowequation.fe_SSA = 'P1bubble'
    md.flowequation.FSNitscheGamma = 1e5
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.isSIA = 0
        self.isL1L2 = 0
        self.isSSA = 0
        self.isMOLHO = 0
        self.isHO = 0
        self.isFS = 0
        self.isNitscheBC = 0
        self.FSNitscheGamma = 1e6
        self.fe_SSA = 'P1'
        self.fe_HO = 'P1'
        self.fe_FS = 'MINIcondensed'
        self.augmented_lagrangian_r = 1
        self.augmented_lagrangian_rhop = 1
        self.augmented_lagrangian_rlambda = 1
        self.augmented_lagrangian_rholambda = 1
        self.XTH_theta = 0
        self.vertex_equation = np.nan
        self.element_equation = np.nan
        self.borderSSA = np.nan
        self.borderHO = np.nan
        self.borderFS = np.nan

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   flow equation parameters:\n'

        s += '{}\n'.format(param_utils.fielddisplay(self, 'isSIA', "is the Shallow Ice Approximation (SIA) used?"))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'isSSA', "is the Shelfy-Stream Approximation (SSA) used?"))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'isL1L2', "are L1L2 equations used?"))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'isMOLHO', "are MOno-layer Higher-Order (MOLHO) equations used?"))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'isHO', "is the Higher-Order (HO) approximation used?"))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'isFS', "are the Full-FS (FS) equations used?"))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'isNitscheBC', "is weakly imposed condition used?"))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'FSNitscheGamma', "Gamma value for the Nitsche term (default: 1e6)"))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'fe_SSA', "Finite Element for SSA: 'P1', 'P1bubble' 'P1bubblecondensed' 'P2'"))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'fe_HO', "Finite Element for HO:  'P1', 'P1bubble', 'P1bubblecondensed', 'P1xP2', 'P2xP1', 'P2', 'P2bubble', 'P1xP3', 'P2xP4'"))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'fe_FS', "Finite Element for FS:  'P1P1' (debugging only) 'P1P1GLS' 'MINIcondensed' 'MINI' 'TaylorHood' 'LATaylorHood' 'XTaylorHood'"))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'vertex_equation', "flow equation for each vertex"))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'element_equation', "flow equation for each element"))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'borderSSA', "vertices on SSA's border (for tiling)"))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'borderHO', "vertices on HO's border (for tiling)"))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'borderFS', "vertices on FS' border (for tiling)"))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - flowequation Class'
        return s

    # Marshall method for saving the flowequation parameters
    def marshall_class(self, fid, prefix, md = None):
        """
        Marshall [flowequation] parameters to a binary file.

        Parameters
        ----------
        fid : file object
            The file object to write the binary data to.
        prefix : str
            Prefix string used for data identification in the binary file.
        md : ISSM model object, optional.
            ISSM model object needed in some cases.

        Returns
        -------
        None
        """
        ## Write Boolean fields
        fieldnames = ['isSIA', 'isSSA', 'isL1L2', 'isMOLHO', 'isHO', 'isFS', 'isNitscheBC']
        for field in fieldnames:
            execute.WriteData(fid, prefix, obj = self, fieldname = field, format = 'Boolean')

        ## Write Double fields
        fieldnames = ['FSNitscheGamma', 'augmented_lagrangian_r', 'augmented_lagrangian_rhop',
                  'augmented_lagrangian_rlambda', 'augmented_lagrangian_rholambda']
        for field in fieldnames:
            execute.WriteData(fid, prefix, obj = self, fieldname = field, format = 'Double')

        ## Write String fields
        fieldnames = ['fe_SSA', 'fe_HO', 'fe_FS']
        for field in fieldnames:
            execute.WriteData(fid, prefix, obj = self, fieldname = field, data = getattr(self, field), format = 'String')
        
        ## Write DoubleMat fields (mattype 1)
        fieldnames = ['borderSSA', 'borderHO', 'borderFS']
        for field in fieldnames:
            execute.WriteData(fid, prefix, obj = self, fieldname = field, format = 'DoubleMat', mattype = 1)

        ## Write other fields
        execute.WriteData(fid, prefix, obj = self, fieldname = 'XTH_theta', data = self.XTH_theta, format = 'Double')
        execute.WriteData(fid, prefix, name = 'md.flowequation.vertex_equation', data = self.vertex_equation, format = 'DoubleMat', mattype = 1)
        execute.WriteData(fid, prefix, name = 'md.flowequation.element_equation', data = self.element_equation, format = 'DoubleMat', mattype = 2)