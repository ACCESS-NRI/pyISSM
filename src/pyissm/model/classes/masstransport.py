import numpy as np
from pyissm.model.classes import class_utils
from pyissm.model.classes import class_registry
from pyissm.model import execute, mesh

@class_registry.register_class
class masstransport(class_registry.manage_state):
    """
    Mass transport solution parameters class for ISSM.

    This class encapsulates parameters for configuring mass transport simulations in the ISSM (Ice Sheet System Model) framework.
    It controls ice thickness evolution, free surface behavior, stabilization methods, and hydrostatic adjustments
    for both grounded and floating ice.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    spcthickness : ndarray, default=nan
        Thickness constraints (NaN means no constraint) [m].
    isfreesurface : int, default=0
        Do we use free surfaces (FS only) or mass conservation.
    min_thickness : float, default=1.0
        Minimum ice thickness allowed [m].
    hydrostatic_adjustment : str, default='Absolute'
        Adjustment of ice shelves surface and bed elevations: 'Incremental' or 'Absolute'.
    stabilization : int, default=1
        Stabilization method: 0=no stabilization, 1=artificial diffusion, 2=streamline upwinding, 3=discontinuous Galerkin, 4=flux corrected transport, 5=streamline upwind Petrov-Galerkin (SUPG).
    vertex_pairing : float, default=nan
        Vertex pairing parameter. Used during consistency checks.
    penalty_factor : float, default=3
        Penalty factor for constraint enforcement.
    requested_outputs : list, default=['default']
        Additional outputs requested.

    Methods
    -------
    __init__(self, other=None)
        Initializes the masstransport parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the masstransport parameters.
    __str__(self)
        Returns a short string identifying the class.
    process_outputs(self, md=None, return_default_outputs=False)
        Process requested outputs, expanding 'default' to appropriate outputs.
    marshall_class(self, fid, prefix, md=None)
        Marshall parameters to a binary file

    Examples
    --------
    md.masstransport = pyissm.model.classes.masstransport()
    md.masstransport.min_thickness = 10.0
    md.masstransport.stabilization = 2
    md.masstransport.hydrostatic_adjustment = 'Incremental'
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.spcthickness = np.nan
        self.isfreesurface = 0
        self.min_thickness = 1.
        self.hydrostatic_adjustment = 'Absolute'
        self.stabilization = 1
        self.vertex_pairing = np.nan
        self.penalty_factor = 3
        self.requested_outputs = ['default']

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   Masstransport solution parameters:\n'

        s += '{}\n'.format(class_utils.fielddisplay(self, 'spcthickness', 'thickness constraints (NaN means no constraint) [m]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'isfreesurface', 'do we use free surfaces (FS only) or mass conservation'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'min_thickness', 'minimum ice thickness allowed [m]'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'hydrostatic_adjustment', 'adjustment of ice shelves surface and bed elevations: ''Incremental'' or ''Absolute'' '))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'stabilization', '0: no stabilization, 1: artificial diffusion, 2: streamline upwinding, 3: discontinuous Galerkin, 4: flux corrected transport, 5: streamline upwind Petrov-Galerkin (SUPG)'))
        s += '{}\n'.format(class_utils.fielddisplay(self, 'requested_outputs', 'additional outputs requested'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - masstransport Class'
        return s
    
    # Extrude to 3D mesh
    def extrude(self, md):
        """
        Extrude masstransport fields to 3D
        """
        self.spcthickness = mesh.project_3d(md, vector = self.spcthickness, type = 'node')
            
        return self
    
    # Check model consistency
    def check_consistency(self, md, solution, analyses):
        # Early return if required analyses and solution are not requested
        if ('MasstransportAnalysis' not in analyses) or (solution == 'TransientSolution' and not md.transient.ismasstransport):
            return md

        class_utils.check_field(md, fieldname = 'masstransport.spcthickness', timeseries = True, allow_inf = False)
        class_utils.check_field(md, fieldname = 'masstransport.isfreesurface', values = [0, 1])
        class_utils.check_field(md, fieldname = 'masstransport.hydrostatic_adjustment', values = ['Absolute', 'Incremental'])
        class_utils.check_field(md, fieldname = 'masstransport.stabilization', values = [0, 1, 2, 3, 4, 5])
        class_utils.check_field(md, fieldname = 'masstransport.min_thickness', gt = 0)
        class_utils.check_field(md, fieldname = 'masstransport.requested_outputs', string_list = True)
        if not np.any(np.isnan(self.vertex_pairing)) and len(self.vertex_pairing) > 0:
            class_utils.check_field(md, fieldname = 'stressbalance.vertex_pairing', gt = 0)

        return md
    
    # Process requested outputs, expanding 'default' to appropriate outputs
    def process_outputs(self,
                        md = None,
                        return_default_outputs = False):
        """
        Process requested outputs, expanding 'default' to appropriate outputs.

        Parameters
        ----------
        md : ISSM model object, optional
            Model object containing mesh information.
        return_default_outputs : bool, default=False
            Whether to also return the list of default outputs.
            
        Returns
        -------
        outputs : list
            List of output strings with 'default' expanded to actual output names.
        default_outputs : list, optional
            Returned only if `return_default_outputs=True`.
        """

        outputs = []

        ## Set default_outputs
        default_outputs = ['Thickness', 'Surface', 'Base']

        ## Loop through all requested outputs
        for item in self.requested_outputs:
            
            ## Process default outputs
            if item == 'default':
                    outputs.extend(default_outputs)

            ## Append other requested outputs (not defaults)
            else:
                outputs.append(item)

        if return_default_outputs:
            return outputs, default_outputs
        return outputs
        
    # Marshall method for saving the masstransport parameters
    def marshall_class(self, fid, prefix, md = None):
        """
        Marshall [masstransport] parameters to a binary file.

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

        ## Write header field
        # NOTE: data types must match the expected types in the ISSM code.
        execute.WriteData(fid, prefix, obj = self, fieldname = 'spcthickness', format = 'DoubleMat', mattype = 1, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'isfreesurface', format = 'Boolean')
        execute.WriteData(fid, prefix, obj = self, fieldname = 'min_thickness', format = 'Double')
        execute.WriteData(fid, prefix, obj = self, fieldname = 'hydrostatic_adjustment', format = 'String')
        execute.WriteData(fid, prefix, obj = self, fieldname = 'stabilization', format = 'Integer')
        execute.WriteData(fid, prefix, obj = self, fieldname = 'vertex_pairing', format = 'DoubleMat', mattype = 3)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'penalty_factor', format = 'Double')
        execute.WriteData(fid, prefix, name = 'md.masstransport.requested_outputs', data = self.process_outputs(md), format = 'StringArray')

