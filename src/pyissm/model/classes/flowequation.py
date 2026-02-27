import numpy as np
import warnings
from pyissm.model.classes import class_utils, class_registry
from pyissm.model import execute, mesh

@class_registry.register_class
class flowequation(class_registry.manage_state):
    """
    Flow equation class for ISSM.

    This class contains parameters for configuring flow equations in the ISSM framework.
    It allows users to select and configure different ice flow approximations such as SIA,
    SSA, Higher-Order, and Full-Stokes, along with their associated numerical methods and
    boundary conditions.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in ``other`` differ from default
        values, they will override the default values.

    Attributes
    ----------
    isSIA : :class:`int`, default=0
        Is the Shallow Ice Approximation (SIA) used?
    isL1L2 : :class:`int`, default=0
        Are L1L2 equations used?
    isSSA : :class:`int`, default=0
        Is the Shelfy-Stream Approximation (SSA) used?
    isMOLHO : :class:`int`, default=0
        Are MOno-layer Higher-Order (MOLHO) equations used?
    isHO : :class:`int`, default=0
        Is the Higher-Order (HO) approximation used?
    isFS : :class:`int`, default=0
        Are the Full-Stokes (FS) equations used?
    isNitscheBC : :class:`int`, default=0
        Is weakly imposed condition used?
    FSNitscheGamma : :class:`float`, default=1e6
        Gamma value for the Nitsche term.
    fe_SSA : :class:`str`, default='P1'
        Finite Element for SSA: 'P1', 'P1bubble', 'P1bubblecondensed', 'P2'.
    fe_HO : :class:`str`, default='P1'
        Finite Element for HO: 'P1', 'P1bubble', 'P1bubblecondensed', 'P2'.
    fe_FS : :class:`str`, default='MINIcondensed'
        Finite Element for FS: 'MINI', 'MINIcondensed', 'TaylorHood', 'XTaylorHood', 'OneLayerP4z', 'CrouzeixRaviart'.
    augmented_lagrangian_r : :class:`float`, default=1
        Augmented Lagrangian parameter r.
    augmented_lagrangian_rhop : :class:`float`, default=1
        Augmented Lagrangian parameter rhop.
    augmented_lagrangian_rlambda : :class:`float`, default=1
        Augmented Lagrangian parameter rlambda.
    augmented_lagrangian_rholambda : :class:`float`, default=1
        Augmented Lagrangian parameter rholambda.
    XTH_theta : :class:`float`, default=0
        XTH theta parameter.
    vertex_equation : :class:`float`, default=np.nan
        Vertex equation parameter.
    element_equation : :class:`float`, default=np.nan
        Element equation parameter.
    borderSSA : :class:`float`, default=np.nan
        Border parameter for SSA.
    borderHO : :class:`float`, default=np.nan
        Border parameter for HO.
    borderFS : :class:`float`, default=np.nan
        Border parameter for FS.

    Examples
    --------
    .. code-block:: python
    
        >>> md.flowequation = pyissm.model.classes.flowequation()
        >>> md.flowequation.isSSA = 1
        >>> md.flowequation.fe_SSA = 'P1bubble'
        >>> md.flowequation.FSNitscheGamma = 1e5
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

        s += '{}\n'.format(class_utils._field_display(self, 'isSIA', "is the Shallow Ice Approximation (SIA) used?"))
        s += '{}\n'.format(class_utils._field_display(self, 'isSSA', "is the Shelfy-Stream Approximation (SSA) used?"))
        s += '{}\n'.format(class_utils._field_display(self, 'isL1L2', "are L1L2 equations used?"))
        s += '{}\n'.format(class_utils._field_display(self, 'isMOLHO', "are MOno-layer Higher-Order (MOLHO) equations used?"))
        s += '{}\n'.format(class_utils._field_display(self, 'isHO', "is the Higher-Order (HO) approximation used?"))
        s += '{}\n'.format(class_utils._field_display(self, 'isFS', "are the Full-FS (FS) equations used?"))
        s += '{}\n'.format(class_utils._field_display(self, 'isNitscheBC', "is weakly imposed condition used?"))
        s += '{}\n'.format(class_utils._field_display(self, 'FSNitscheGamma', "Gamma value for the Nitsche term (default: 1e6)"))
        s += '{}\n'.format(class_utils._field_display(self, 'fe_SSA', "Finite Element for SSA: 'P1', 'P1bubble' 'P1bubblecondensed' 'P2'"))
        s += '{}\n'.format(class_utils._field_display(self, 'fe_HO', "Finite Element for HO:  'P1', 'P1bubble', 'P1bubblecondensed', 'P1xP2', 'P2xP1', 'P2', 'P2bubble', 'P1xP3', 'P2xP4'"))
        s += '{}\n'.format(class_utils._field_display(self, 'fe_FS', "Finite Element for FS:  'P1P1' (debugging only) 'P1P1GLS' 'MINIcondensed' 'MINI' 'TaylorHood' 'LATaylorHood' 'XTaylorHood'"))
        s += '{}\n'.format(class_utils._field_display(self, 'vertex_equation', "flow equation for each vertex"))
        s += '{}\n'.format(class_utils._field_display(self, 'element_equation', "flow equation for each element"))
        s += '{}\n'.format(class_utils._field_display(self, 'borderSSA', "vertices on SSA's border (for tiling)"))
        s += '{}\n'.format(class_utils._field_display(self, 'borderHO', "vertices on HO's border (for tiling)"))
        s += '{}\n'.format(class_utils._field_display(self, 'borderFS', "vertices on FS' border (for tiling)"))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - flowequation Class'
        return s
        
    # Extrude to 3D mesh
    def _extrude(self, md):
        """
        Extrude flowequation fields to 3D
        """
        self.element_equation = mesh._project_3d(md, vector = self.element_equation, type = 'element')
        self.vertex_equation = mesh._project_3d(md, vector = self.vertex_equation, type = 'node')
        self.borderSSA = mesh._project_3d(md, vector = self.borderSSA, type = 'node')
        self.borderHO = mesh._project_3d(md, vector = self.borderHO, type = 'node')
        self.borderFS = mesh._project_3d(md, vector = self.borderFS, type = 'node')
            
        return self
    
    # Check model consistency
    def check_consistency(self, md, solution, analyses):
        """
        Check consistency of the [flowequation.default] parameters.

        Parameters
        ----------
        md : :class:`pyissm.model.Model`
            The model object to check.
        solution : :class:`pyissm.model.solution`
            The solution object to check.
        analyses : list of :class:`str`
            List of analyses to check consistency for.

        Returns 
        -------
        md : :class:`pyissm.model.Model`
            The model object with any consistency errors noted.
        """

        # Early return if necessary analyses and solutions not present
        if ('StressbalanceAnalysis' not in analyses and 'StressbalanceSIAAnalysis' not in analyses) or (solution == 'TransientSolution' and not md.transient.isstressbalance):
            return md

        class_utils.check_field(md, fieldname = "flowequation.isSIA", scalar = True, values = [0, 1])
        class_utils.check_field(md, fieldname = "flowequation.isSSA", scalar = True, values = [0, 1])
        class_utils.check_field(md, fieldname = "flowequation.isL1L2", scalar = True, values = [0, 1])
        class_utils.check_field(md, fieldname = "flowequation.isMOLHO", scalar = True, values = [0, 1])
        class_utils.check_field(md, fieldname = "flowequation.isHO", scalar = True, values = [0, 1])
        class_utils.check_field(md, fieldname = "flowequation.isFS", scalar = True, values = [0, 1])
        class_utils.check_field(md, fieldname = "flowequation.isNitscheBC", scalar = True, values = [0, 1])
        class_utils.check_field(md, fieldname = "flowequation.FSNitscheGamma", scalar = True, ge = 0.0)
        class_utils.check_field(md, fieldname = 'flowequation.fe_SSA', values = ['P1', 'P1bubble', 'P1bubblecondensed', 'P2', 'P2bubble'])
        class_utils.check_field(md, fieldname = 'flowequation.fe_HO', values =  ['P1', 'P1bubble', 'P1bubblecondensed', 'P1xP2', 'P2xP1', 'P2', 'P2bubble', 'P1xP3', 'P2xP4'])
        class_utils.check_field(md, fieldname = 'flowequation.fe_FS', values = ['P1P1', 'P1P1GLS', 'MINIcondensed', 'MINI', 'TaylorHood', 'LATaylorHood', 'XTaylorHood', 'OneLayerP4z', 'CrouzeixRaviart', 'LACrouzeixRaviart'])
        class_utils.check_field(md, fieldname = 'flowequation.borderSSA', size = (md.mesh.numberofvertices, ), values = [0, 1])
        class_utils.check_field(md, fieldname = 'flowequation.borderHO', size = (md.mesh.numberofvertices, ), values = [0, 1])
        class_utils.check_field(md, fieldname = 'flowequation.borderFS', size = (md.mesh.numberofvertices, ), values = [0, 1])
        class_utils.check_field(md, fieldname = "flowequation.augmented_lagrangian_r", scalar = True, gt = 0.0)
        class_utils.check_field(md, fieldname = "flowequation.augmented_lagrangian_rhop", scalar = True, gt = 0.0)
        class_utils.check_field(md, fieldname = "flowequation.augmented_lagrangian_rlambda", scalar = True, gt = 0.0)
        class_utils.check_field(md, fieldname = "flowequation.augmented_lagrangian_rholambda", scalar = True, gt = 0.0)
        class_utils.check_field(md, fieldname = "flowequation.XTH_theta", scalar = True, ge = 0.0, lt = 0.5)

        if md.mesh.domain_type() == '2Dhorizontal':
            class_utils.check_field(md, fieldname = 'flowequation.vertex_equation', size = (md.mesh.numberofvertices, ), values = [1, 2, 4])
            class_utils.check_field(md, fieldname = 'flowequation.element_equation', size = (md.mesh.numberofelements, ), values = [1, 2, 4])
        elif md.mesh.domain_type() == '3Dsurface':
            class_utils.check_field(md, fieldname = 'flowequation.vertex_equation', size = (md.mesh.numberofvertices, ), values = np.arange(1, 2 + 1))
            class_utils.check_field(md, fieldname = 'flowequation.element_equation', size = (md.mesh.numberofelements, ), values = np.arange(1, 2 + 1))
        elif md.mesh.domain_type() == '2Dvertical':
            class_utils.check_field(md, fieldname = 'flowequation.vertex_equation', size = (md.mesh.numberofvertices, ), values = [2, 5, 6])
            class_utils.check_field(md, fieldname = 'flowequation.element_equation', size = (md.mesh.numberofelements, ), values = [2, 5, 6])
        elif md.mesh.domain_type() == '3D':
            class_utils.check_field(md, fieldname = 'flowequation.vertex_equation', size = (md.mesh.numberofvertices, ), values = np.arange(0, 9 + 1))
            class_utils.check_field(md, fieldname = 'flowequation.element_equation', size = (md.mesh.numberofelements, ), values = np.arange(0, 9 + 1))
        else:
            raise RuntimeError(f'pyissm.model.classes.flowequation.check_consistency: unknown domain_type: {md.mesh.domain_type()}.')
    
        if not (self.isSIA or self.isSSA or self.isL1L2 or self.isMOLHO or self.isHO or self.isFS):
            md.check_message("no element types set for this model")
        if 'StressbalanceSIAAnalysis' in analyses:
            if any(self.element_equation == 1):
                if np.any(np.logical_and(self.vertex_equation, md.mask.ocean_levelset)):
                    warnings.warn("\n !!! Warning: SIA's model is not consistent on ice shelves !!!\n")
        
        return md

    # Marshall method for saving the flowequation parameters
    def marshall_class(self, fid, prefix, md  =  None):
        """
        Marshall [flowequation.default] parameters to a binary file.

        Parameters
        ----------
        fid : :class:`file object`
            The file object to write the binary data to.
        prefix : :class:`str`
            Prefix string used for data identification in the binary file.
        md : :class:`pyissm.model.Model`, optional
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