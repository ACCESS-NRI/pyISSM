import numpy as np
from pyissm.model.classes import class_utils, class_registry
from pyissm.model import execute

@class_registry.register_class
class esa(class_registry.manage_state):
    """
    Elastic solid earth adjustment (ESA) class for ISSM.

    This class contains parameters for elastic solid earth adjustment in the ISSM framework.
    ESA models the instantaneous elastic response of the solid Earth to ice loading changes,
    which is important for studying present-day ice sheet mass balance and sea level change.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in ``other`` differ from default
        values, they will override the default values.

    Attributes
    ----------
    deltathickness : :class:`float`, default=np.nan
        Thickness change: ice height equivalent [m].
    love_h : :class:`float`, default=0.0
        Load Love number for radial displacement.
    love_l : :class:`float`, default=0.0
        Load Love number for horizontal displacements.
    hemisphere : :class:`float`, default=0.0
        North-south, East-west components of 2-D horizontal displacement vector: -1 south, 1 north.
    degacc : :class:`float`, default=0.01
        Accuracy (default 0.01 deg) for numerical discretization of the Green's functions.
    requested_outputs : :class:`list`, default=['default']
        Additional outputs requested (default: EsaUmotion).
    transitions : :class:`list`, default=[]
        Indices into parts of the mesh that will be icecaps.

    Examples
    --------
    .. code-block:: python
    
        >>> md.esa = pyissm.model.classes.esa()
        >>> md.esa.deltathickness = thickness_change
        >>> md.esa.love_h = 0.6
        >>> md.esa.love_l = 0.1
        >>> md.esa.degacc = 0.005
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.deltathickness = np.nan
        self.love_h = 0.
        self.love_l = 0.
        self.hemisphere = 0.
        self.degacc = 0.01
        self.requested_outputs = ['default']
        self.transitions = []

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   esa parameters:\n'

        s += '{}\n'.format(class_utils._field_display(self, 'deltathickness', 'thickness change: ice height equivalent [m]'))
        s += '{}\n'.format(class_utils._field_display(self, 'love_h', 'load Love number for radial displacement'))
        s += '{}\n'.format(class_utils._field_display(self, 'love_l', 'load Love number for horizontal displaements'))
        s += '{}\n'.format(class_utils._field_display(self, 'hemisphere', 'North-south, East-west components of 2-D horiz displacement vector:-1 south, 1 north'))
        s += '{}\n'.format(class_utils._field_display(self, 'degacc', 'accuracy (default .01 deg) for numerical discretization of the Green''s functions'))
        s += '{}\n'.format(class_utils._field_display(self, 'transitions', 'indices into parts of the mesh that will be icecaps'))
        s += '{}\n'.format(class_utils._field_display(self, 'requested_outputs', 'additional outputs requested (default: EsaUmotion)'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - esa Class'
        return s
    
    # Check model consistency
    def check_consistency(self, md, solution, analyses):
        """
        Check consistency of the [esa] parameters.

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

        # Early return if EsaAnalysis not specified
        if (solution != 'EsaAnalysis'):
            return md
        
        class_utils.check_field(md, fieldname = "esa.deltathickness", allow_nan = True, allow_inf = True, size = (md.mesh.numberofelements, ))
        class_utils.check_field(md, fieldname = "esa.love_h", allow_nan = False, allow_inf = False)
        class_utils.check_field(md, fieldname = "esa.love_l", allow_nan = False, allow_inf = False)
        class_utils.check_field(md, fieldname = "esa.hemisphere", allow_nan = False, allow_inf = False)
        class_utils.check_field(md, fieldname = "esa.degacc", size = (1, 1), ge = 1e-10)
        class_utils.check_field(md, fieldname = "esa.requested_outputs", string_list = True)

        if (np.size(self.love_h, 1) != np.size(self.love_l, 0)):
            raise ValueError('pyissm.model.classes.esa.check_consistency: love_h and love_l must be the same size.')


    # Process requested outputs, expanding 'default' to appropriate outputs
    def process_outputs(self,
                        md = None,
                        return_default_outputs = False):
        """
        Process requested outputs for [esa] parameters, expanding 'default' to appropriate outputs.

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
        default_outputs = ['EsaUmotion']

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

    # Marshall method for saving the esa parameters
    def marshall_class(self, fid, prefix, md = None):
        """
        Marshall [esa] parameters to a binary file.

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
                
        ## Write fields
        execute.WriteData(fid, prefix, obj = self, fieldname = 'deltathickness', format = 'DoubleMat', mattype = 2)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'love_h', format = 'DoubleMat', mattype = 1)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'love_l', format = 'DoubleMat', mattype = 1)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'hemisphere', format = 'Integer')
        execute.WriteData(fid, prefix, obj = self, fieldname = 'degacc', format = 'Double')
        execute.WriteData(fid, prefix, obj = self, fieldname = 'transitions', format = 'MatArray')
        execute.WriteData(fid, prefix, name = 'md.esa.requested_outputs', data = self.process_outputs(md), format = 'StringArray')