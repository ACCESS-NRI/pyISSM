import numpy as np
from . import param_utils
from . import class_registry
from .. import execute

@class_registry.register_class
class groundingline(class_registry.manage_state):
    """
    Grounding line migration parameters class for ISSM.

    This class encapsulates parameters for configuring grounding line migration in the ISSM (Ice Sheet System Model) framework.
    It controls how the grounding line (boundary between grounded and floating ice) moves during simulations,
    including migration methods and interpolation schemes for friction and melting on partially floating elements.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    migration : str, default='SubelementMigration'
        Type of grounding line migration: 'SoftMigration', 'SubelementMigration', 'AggressiveMigration', 'Contact', 'None'.
    friction_interpolation : str, default='SubelementFriction1'
        Type of friction interpolation on partially floating elements: 'SubelementFriction1', 'SubelementFriction2', 'NoFrictionOnPartiallyFloating'.
    melt_interpolation : str, default='NoMeltOnPartiallyFloating'
        Type of melt interpolation on partially floating elements: 'SubelementMelt1', 'SubelementMelt2', 'IntrusionMelt', 'NoMeltOnPartiallyFloating', 'FullMeltOnPartiallyFloating'.
    intrusion_distance : float, default=0
        Distance of seawater intrusion from grounding line [m].
    requested_outputs : list, default=['default']
        Additional outputs requested for grounding line analysis.

    Methods
    -------
    __init__(self, other=None)
        Initializes the groundingline parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the groundingline parameters.
    __str__(self)
        Returns a short string identifying the class.
    process_outputs(self, md=None, return_default_outputs=False)
        Process requested outputs, expanding 'default' to appropriate outputs.
    marshall_class(self, fid, prefix, md=None)
        Marshall parameters to a binary file.

    Examples
    --------
    md.groundingline = pyissm.param.groundingline()
    md.groundingline.migration = 'AggressiveMigration'
    md.groundingline.friction_interpolation = 'SubelementFriction2'
    md.groundingline.melt_interpolation = 'SubelementMelt1'
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.migration = 'SubelementMigration'
        self.friction_interpolation = 'SubelementFriction1'
        self.melt_interpolation = 'NoMeltOnPartiallyFloating'
        self.intrusion_distance = 0
        self.requested_outputs = ['default']

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   grounding line migration parameters:\n'

        s += '{}\n'.format(param_utils.fielddisplay(self, 'migration', 'type of grounding line migration: \'SoftMigration\', \'SubelementMigration\', \'AggressiveMigration\', \'Contact\', \'None\''))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'friction_interpolation', 'type of friction interpolation on partially floating elements: ''SubelementFriction1'', ''SubelementFriction2'', ''NoFrictionOnPartiallyFloating'''))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'melt_interpolation', 'type of melt interpolation on partially floating elements: \'SubelementMelt1\', \'SubelementMelt2\', \'IntrusionMelt\', \'NoMeltOnPartiallyFloating\', \'FullMeltOnPartiallyFloating\''))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'intrusion_distance', 'distance of seawater intrusion from grounding line [m]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'requested_outputs', 'additional outputs requested'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - groundingline Class'
        return s
    
    # Check model consistency
    def check_consistency(self, md, solution, analyses):  # {{{
        param_utils.check_field(md, fieldname = 'groundingline.migration', values = ['None', 'SubelementMigration', 'AggressiveMigration', 'SoftMigration', 'Contact', 'GroundingOnly'])
        param_utils.check_field(md, fieldname = 'groundingline.friction_interpolation', values = ['SubelementFriction1', 'SubelementFriction2', 'NoFrictionOnPartiallyFloating'])
        param_utils.check_field(md, fieldname = 'groundingline.melt_interpolation', values = ['SubelementMelt1', 'SubelementMelt2', 'IntrusionMelt', 'NoMeltOnPartiallyFloating', 'FullMeltOnPartiallyFloating'])
        param_utils.check_field(md, fieldname = 'groundingline.intrusion_distance', ge = 0)
        param_utils.check_field(md, fieldname = 'groundingline.requested_outputs', string_list = True)

        if(not self.migration == 'None' and md.transient.isgroundingline and solution == 'TransientSolution'):
            if np.any(np.isnan(md.geometry.bed)):
                md.checkmessage("requesting grounding line migration, but bathymetry is absent!")
            pos = np.nonzero(md.mask.ocean_levelset > 0.)[0]
            if any(np.abs(md.geometry.base[pos] - md.geometry.bed[pos]) > pow(10, -10)):
                md.checkmessage("base not equal to bed on grounded ice!")
            if any(md.geometry.bed - md.geometry.base > pow(10, -9)):
                md.checkmessage("bed superior to base on floating ice!")
                
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
        default_outputs = ['Surface', 'Base','MaskOceanLevelset']

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

    # Marshall method for saving the groundingline parameters
    def marshall_class(self, fid, prefix, md = None):
        """
        Marshall [groundingline] parameters to a binary file.

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

        ## Write String fields
        fieldnames = ['migration', 'friction_interpolation', 'melt_interpolation']
        for fieldname in fieldnames:
            execute.WriteData(fid, prefix, name = 'md.groundingline.' + fieldname, data = getattr(self, fieldname), format = 'String')
        
        ## Write other fields
        execute.WriteData(fid, prefix, obj = self, fieldname = 'intrusion_distance', format = 'Double')
        execute.WriteData(fid, prefix, name = 'md.groundingline.requested_outputs', data = self.process_outputs(md), format = 'StringArray')