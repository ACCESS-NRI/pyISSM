import numpy as np
from pyissm.model.classes import class_utils, class_registry
from pyissm.model import execute

@class_registry.register_class
class groundingline(class_registry.manage_state):
    """
    Grounding line migration class for ISSM.

    This class contains parameters for configuring grounding line migration in the ISSM framework.
    It controls how the grounding line (boundary between grounded and floating ice) moves during simulations,
    including migration methods and interpolation schemes for friction and melting on partially floating elements.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in ``other`` differ from default
        values, they will override the default values.

    Attributes
    ----------
    migration : :class:`str`, default='SubelementMigration'
        Type of grounding line migration: 'SoftMigration', 'SubelementMigration', 'AggressiveMigration', 'Contact', 'None'.
    friction_interpolation : :class:`str`, default='SubelementFriction1'
        Type of friction interpolation on partially floating elements: 'SubelementFriction1', 'SubelementFriction2', 'NoFrictionOnPartiallyFloating'.
    melt_interpolation : :class:`str`, default='NoMeltOnPartiallyFloating'
        Type of melt interpolation on partially floating elements: 'SubelementMelt1', 'SubelementMelt2', 'IntrusionMelt', 'NoMeltOnPartiallyFloating', 'FullMeltOnPartiallyFloating'.
    intrusion_distance : :class:`float`, default=0
        Distance of seawater intrusion from grounding line [m].
    requested_outputs : :class:`list`, default=['default']
        Additional outputs requested for grounding line analysis.

    Examples
    --------
    .. code-block:: python

        >>> md.groundingline = pyissm.model.classes.groundingline()
        >>> md.groundingline.migration = 'AggressiveMigration'
        >>> md.groundingline.friction_interpolation = 'SubelementFriction2'
        >>> md.groundingline.melt_interpolation = 'SubelementMelt1'
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

        s += '{}\n'.format(class_utils._field_display(self, 'migration', 'type of grounding line migration: \'SoftMigration\', \'SubelementMigration\', \'AggressiveMigration\', \'Contact\', \'None\''))
        s += '{}\n'.format(class_utils._field_display(self, 'friction_interpolation', 'type of friction interpolation on partially floating elements: ''SubelementFriction1'', ''SubelementFriction2'', ''NoFrictionOnPartiallyFloating'''))
        s += '{}\n'.format(class_utils._field_display(self, 'melt_interpolation', 'type of melt interpolation on partially floating elements: \'SubelementMelt1\', \'SubelementMelt2\', \'IntrusionMelt\', \'NoMeltOnPartiallyFloating\', \'FullMeltOnPartiallyFloating\''))
        s += '{}\n'.format(class_utils._field_display(self, 'intrusion_distance', 'distance of seawater intrusion from grounding line [m]'))
        s += '{}\n'.format(class_utils._field_display(self, 'requested_outputs', 'additional outputs requested'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - groundingline Class'
        return s
    
    # Check model consistency
    def check_consistency(self, md, solution, analyses):
        """
        Check consistency of the [groundingline] parameters.

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

        class_utils._check_field(md, fieldname = 'groundingline.migration', values = ['None', 'SubelementMigration', 'AggressiveMigration', 'SoftMigration', 'Contact', 'GroundingOnly'])
        class_utils._check_field(md, fieldname = 'groundingline.friction_interpolation', values = ['SubelementFriction1', 'SubelementFriction2', 'NoFrictionOnPartiallyFloating'])
        class_utils._check_field(md, fieldname = 'groundingline.melt_interpolation', values = ['NoMeltOnPartiallyFloating', 'FullMeltOnPartiallyFloating', 'SubelementMelt1', 'SubelementMelt2', 'IntrusionMelt'])
        class_utils._check_field(md, fieldname = 'groundingline.intrusion_distance', ge = 0, allow_nan = False, allow_inf = False)
        class_utils._check_field(md, fieldname = 'groundingline.requested_outputs', string_list = True)

        if(not self.migration == 'None' and md.transient.isgroundingline and solution == 'TransientSolution'):
            if np.any(np.isnan(md.geometry.bed)):
                md.check_message("requesting grounding line migration, but bathymetry is absent!")
            pos = np.nonzero(md.mask.ocean_levelset > 0.)[0]
            if any(np.abs(md.geometry.base[pos] - md.geometry.bed[pos]) > pow(10, -10)):
                md.check_message("base not equal to bed on grounded ice!")
            if any(md.geometry.bed - md.geometry.base > pow(10, -9)):
                md.check_message("bed superior to base on floating ice!")
                
        return md
    
    # Process requested outputs, expanding 'default' to appropriate outputs
    def process_outputs(self,
                        md = None,
                        return_default_outputs = False):
        """
        Process requested outputs for [groundingline] parameters, expanding 'default' to appropriate outputs.

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

        ## Write String fields
        fieldnames = ['migration', 'friction_interpolation', 'melt_interpolation']
        for fieldname in fieldnames:
            execute._write_model_field(fid, prefix, name = 'md.groundingline.' + fieldname, data = getattr(self, fieldname), format = 'String')
        
        ## Write other fields
        execute._write_model_field(fid, prefix, obj = self, fieldname = 'intrusion_distance', format = 'DoubleMat', mattype = 1)
        execute._write_model_field(fid, prefix, name = 'md.groundingline.requested_outputs', data = self.process_outputs(md), format = 'StringArray')