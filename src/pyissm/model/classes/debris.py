import numpy as np
from pyissm.model.classes import class_utils, class_registry
from pyissm.model import execute, mesh

@class_registry.register_class
class debris(class_registry.manage_state):
    """
    Debris transport class for ISSM.

    This class contains parameters for debris transport modeling in the ISSM framework.
    Debris transport simulates the movement and accumulation of rock debris on glacier surfaces,
    which affects surface albedo, melting rates, and overall ice dynamics.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in ``other`` differ from default
        values, they will override the default values.

    Attributes
    ----------
    spcthickness : :class:`numpy.ndarray`, default=np.nan
        Debris thickness constraints (NaN means no constraint) [m].
    min_thickness : :class:`float`, default=0.0
        Minimum debris thickness allowed [m].
    stabilization : :class:`int`, default=2
        Stabilization method: 0=no stabilization, 1=artificial diffusion, 2=streamline upwinding, 3=streamline upwind Petrov-Galerkin (SUPG).
    packingfraction : :class:`float`, default=0.01
        Fraction of debris covered in the ice.
    removalmodel : :class:`int`, default=0
        Frontal removal of debris: 0=no removal, 1=Slope-triggered debris removal, 2=driving-stress triggered debris removal.
    displacementmodel : :class:`int`, default=0
        Debris displacement: 0=no displacement, 1=...
    max_displacementvelocity : :class:`float`, default=0.0
        Maximum velocity of debris transport (v_ice + v_displacement) [m/a].
    removal_slope_threshold : :class:`float`, default=0.0
        Critical slope [degrees] for removalmodel (1).
    removal_stress_threshold : :class:`float`, default=0.0
        Critical stress [Pa] for removalmodel (2).
    vertex_pairing : :class:`float`, default=np.nan
        Pairs of vertices that are penalized.
    requested_outputs : :class:`list`, default=['default']
        Additional outputs requested.

    Examples
    --------
    .. code-block:: python
    
        >>> md.debris = pyissm.model.classes.debris()
        >>> md.debris.min_thickness = 0.001
        >>> md.debris.packingfraction = 0.02
        >>> md.debris.stabilization = 2
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

        s += '{}\n'.format(class_utils._field_display(self,'spcthickness','debris thickness constraints (NaN means no constraint) [m]'))
        s += '{}\n'.format(class_utils._field_display(self,'min_thickness','minimum debris thickness allowed [m]'))
        s += '{}\n'.format(class_utils._field_display(self,'packingfraction','fraction of debris covered in the ice'))
        s += '{}\n'.format(class_utils._field_display(self,'stabilization','0: no stabilization, 1: artificial diffusion, 2: streamline upwinding, 3: streamline upwind Petrov-Galerkin (SUPG)'))
        s += '{}\n'.format(class_utils._field_display(self,'removalmodel','frontal removal of debris. 0: no removal, 1: Slope-triggered debris removal, 2: driving-stress triggered debris removal'))
        s += '{}\n'.format(class_utils._field_display(self,'displacementmodel','debris displacement. 0: no displacement, 1: ...'))
        s += '{}\n'.format(class_utils._field_display(self,'max_displacementvelocity','maximum velocity of debris transport (v_ice + v_displacement) (m/a)'))
        s += '{}\n'.format(class_utils._field_display(self,'removal_slope_threshold','critical slope (degrees) for removalmodel (1)'))
        s += '{}\n'.format(class_utils._field_display(self,'removal_stress_threshold','critical stress (Pa) for removalmodel (2)'))

        s += '\n      {}\n'.format('Penalty options:')
        s += '{}\n'.format(class_utils._field_display(self,'vertex_pairing','pairs of vertices that are penalized'))
        s += '{}\n'.format(class_utils._field_display(self,'requested_outputs','additional outputs requested'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - debris Class'
        return s
    
    # Extrude to 3D mesh
    def _extrude(self, md):
        """
        Extrude debris fields to 3D
        """
        self.spcthickness = mesh._project_3d(md, vector = self.spcthickness, type = 'node')
            
        return self

    # Check model consistency
    def check_consistency(self, md, solution, analyses):
        """
        Check consistency of the [debris] parameters.

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

        # Early return if not mass transport analysis or transient with debris transport
        if not 'MassTransportAnalysis' in analyses or solution == 'TransientSolution' and not md.transient.isdebris:
            return md

        class_utils.check_field(md, fieldname = "debris.spcthickness")
        class_utils.check_field(md, fieldname = "debris.stabilization", scalar = True, values = [0, 1, 2, 3, 4, 5])
        class_utils.check_field(md, fieldname = "debris.min_thickness", ge = 0)
        class_utils.check_field(md, fieldname = "debris.packingfraction", ge = 0)
        class_utils.check_field(md, fieldname = "debris.removalmodel", values = [0, 1, 2])
        class_utils.check_field(md, fieldname = "debris.displacementmodel", values = [0, 1, 2])
        class_utils.check_field(md, fieldname = "debris.max_displacementvelocity", ge = 0)
        class_utils.check_field(md, fieldname = "debris.removal_slope_threshold", ge = 0)
        class_utils.check_field(md, fieldname = "debris.removal_stress_threshold", ge = 0)
        class_utils.check_field(md, fieldname = "debris.requested_outputs", string_list = True)
        if not np.any(np.isnan(md.stressbalance.vertex_pairing)) and len(md.stressbalance.vertex_pairing) > 0:
            md = class_utils.check_field(md, fieldname = "stressbalance.vertex_pairing", gt = 0)
            
        return md
    

    # Process requested outputs, expanding 'default' to appropriate outputs
    def process_outputs(self,
                        md = None,
                        return_default_outputs = False):
        """
        Process requested outputs for [debris] parameters, expanding 'default' to appropriate outputs.

        Parameters
        ----------
        md : :class:`pyissm.model.Model`, optional
            Model object containing mesh information.
        return_default_outputs : :class:`bool`, default=False
            Whether to also return the list of default outputs.
            
        Returns
        -------
        outputs : :class:`list`
            List of output strings with 'default' expanded to actual output names.
        default_outputs : :class:`list`, optional
            Returned only if `return_default_outputs=True`.
        """

        outputs = []

        ## Set default_outputs
        default_outputs = ['DebrisThickness', 'DebrisMaskNodeActivation', 'VxDebris', 'VyDebris']

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
    
    # Marshall method for saving the debris parameters
    def marshall_class(self, fid, prefix, md = None):
        """
        Marshall [debris] parameters to a binary file.

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
        
        ## Write Integer fields
        fieldnames = ['stabilization', 'removalmodel', 'displacementmodel']
        for field in fieldnames:
            execute._write_model_field(fid, prefix, obj = self, fieldname = field, format = 'Integer')
        
        ## Write Double fields
        fieldnames = ['max_displacementvelocity', 'removal_slope_threshold', 'removal_stress_threshold', 'packingfraction', 'min_thickness']
        for field in fieldnames:
            execute._write_model_field(fid, prefix, obj = self, fieldname = field, format = 'Double')
    
        ## Write other fields
        execute._write_model_field(fid, prefix, obj = self, fieldname = 'spcthickness', format = 'DoubleMat', mattype = 1, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
        execute._write_model_field(fid, prefix, obj = self, fieldname = 'vertex_pairing', format = 'DoubleMat', mattype = 3)
        execute._write_model_field(fid, prefix, name = 'md.debris.requested_outputs', data = self.process_outputs(md), format = 'StringArray')

