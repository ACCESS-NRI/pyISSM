import numpy as np
from . import param_utils
from . import class_registry
from .. import execute

@class_registry.register_class
class debris(class_registry.manage_state):
    """
    Debris transport parameters class for ISSM.

    This class encapsulates parameters for debris transport modeling in the ISSM (Ice Sheet System Model) framework.
    Debris transport simulates the movement and accumulation of rock debris on glacier surfaces,
    which affects surface albedo, melting rates, and overall ice dynamics.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    spcthickness : ndarray, default=nan
        Debris thickness constraints (NaN means no constraint) [m].
    min_thickness : float, default=0.0
        Minimum debris thickness allowed [m].
    stabilization : int, default=2
        Stabilization method: 0=no stabilization, 1=artificial diffusion, 2=streamline upwinding, 3=streamline upwind Petrov-Galerkin (SUPG).
    packingfraction : float, default=0.01
        Fraction of debris covered in the ice.
    removalmodel : int, default=0
        Frontal removal of debris: 0=no removal, 1=Slope-triggered debris removal, 2=driving-stress triggered debris removal.
    displacementmodel : int, default=0
        Debris displacement: 0=no displacement, 1=...
    max_displacementvelocity : float, default=0.0
        Maximum velocity of debris transport (v_ice + v_displacement) [m/a].
    removal_slope_threshold : float, default=0.0
        Critical slope [degrees] for removalmodel (1).
    removal_stress_threshold : float, default=0.0
        Critical stress [Pa] for removalmodel (2).
    vertex_pairing : float, default=nan
        Pairs of vertices that are penalized.
    requested_outputs : str, default='List of requested outputs'
        Additional outputs requested.

    Methods
    -------
    __init__(self, other=None)
        Initializes the debris parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the debris parameters.
    __str__(self)
        Returns a short string identifying the class.
    process_outputs(self, md)
        Process requested outputs, expanding 'default' to appropriate outputs.
    marshall_class(self, fid, prefix, md=None)
        Marshall parameters to a binary file.

    Examples
    --------
    md.debris = pyissm.param.debris()
    md.debris.min_thickness = 0.001
    md.debris.packingfraction = 0.02
    md.debris.stabilization = 2
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.spcthickness = np.nan
        self.min_thickness = 0.
        self.stabilization = 2
        self.packingfraction = 0.01
        self.removalmodel = 0
        self.displacementmodel = 0
        self.max_displacementvelocity = 0.
        self.removal_slope_threshold = 0.
        self.removal_stress_threshold = 0.
        self.vertex_pairing = np.nan
        self.requested_outputs = ['default']

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   debris solution parameters:\n'

        s += '{}\n'.format(param_utils.fielddisplay(self,'spcthickness','debris thickness constraints (NaN means no constraint) [m]'))
        s += '{}\n'.format(param_utils.fielddisplay(self,'min_thickness','minimum debris thickness allowed [m]'))
        s += '{}\n'.format(param_utils.fielddisplay(self,'packingfraction','fraction of debris covered in the ice'))
        s += '{}\n'.format(param_utils.fielddisplay(self,'stabilization','0: no stabilization, 1: artificial diffusion, 2: streamline upwinding, 3: streamline upwind Petrov-Galerkin (SUPG)'))
        s += '{}\n'.format(param_utils.fielddisplay(self,'removalmodel','frontal removal of debris. 0: no removal, 1: Slope-triggered debris removal, 2: driving-stress triggered debris removal'))
        s += '{}\n'.format(param_utils.fielddisplay(self,'displacementmodel','debris displacement. 0: no displacement, 1: ...'))
        s += '{}\n'.format(param_utils.fielddisplay(self,'max_displacementvelocity','maximum velocity of debris transport (v_ice + v_displacement) (m/a)'))
        s += '{}\n'.format(param_utils.fielddisplay(self,'removal_slope_threshold','critical slope (degrees) for removalmodel (1)'))
        s += '{}\n'.format(param_utils.fielddisplay(self,'removal_stress_threshold','critical stress (Pa) for removalmodel (2)'))

        s += '\n      {}\n'.format('Penalty options:')
        s += '{}\n'.format(param_utils.fielddisplay(self,'vertex_pairing','pairs of vertices that are penalized'))
        s += '{}\n'.format(param_utils.fielddisplay(self,'requested_outputs','additional outputs requested'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - debris Class'
        return s

# Process requested outputs, expanding 'default' to appropriate outputs
    def process_outputs(self, md = None):
        """
        Process requested outputs, expanding 'default' to appropriate outputs.

        Parameters
        ----------
        md : ISSM model object, optional
            Model object containing mesh information.
            
        Returns
        -------
        outputs
            List of output strings with 'default' expanded to actual output names.
        """
        outputs = []
        default_outputs = ['DebrisThickness', 'DebrisMaskNodeActivation', 'VxDebris', 'VyDebris']

        ## Loop through all requested outputs
        for item in self.requested_outputs:
            
            ## Process default outputs
            if item == 'default':
                    outputs.extend(default_outputs)

            ## Append other requested outputs (not defaults)
            else:
                outputs.append(item)

        return outputs
    
    # Marshall method for saving the debris parameters
    def marshall_class(self, fid, prefix, md = None):
        """
        Marshall [debris] parameters to a binary file.

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
        
        ## Write Integer fields
        fieldnames = ['stabilization', 'removalmodel', 'displacementmodel']
        for field in fieldnames:
            execute.WriteData(fid, prefix, obj = self, fieldname = field, format = 'Integer')
        
        ## Write Double fields
        fieldnames = ['max_displacementvelocity', 'removal_slope_threshold', 'removal_stress_threshold', 'packingfraction']
        for field in fieldnames:
            execute.WriteData(fid, prefix, obj = self, fieldname = field, format = 'Double')
    
        ## Write other fields
        execute.WriteData(fid, prefix, obj = self, fieldname = 'spcthickness', format = 'DoubleMat', mattype = 1, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'vertex_pairing', format = 'DoubleMat', mattype = 3)
        execute.WriteData(fid, prefix, name = 'md.debris.requested_outputs', data = self.process_outputs(md), format = 'StringArray')

