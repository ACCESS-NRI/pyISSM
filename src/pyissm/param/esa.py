import numpy as np
from . import param_utils
from . import class_registry
from .. import execute

@class_registry.register_class
class esa(class_registry.manage_state):
    """
    Elastic solid earth adjustment (ESA) parameters class for ISSM.

    This class encapsulates parameters for elastic solid earth adjustment in the ISSM (Ice Sheet System Model) framework.
    ESA models the instantaneous elastic response of the solid Earth to ice loading changes,
    which is important for studying present-day ice sheet mass balance and sea level change.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    deltathickness : float, default=nan
        Thickness change: ice height equivalent [m].
    love_h : float, default=0.0
        Load Love number for radial displacement.
    love_l : float, default=0.0
        Load Love number for horizontal displacements.
    hemisphere : float, default=0.0
        North-south, East-west components of 2-D horizontal displacement vector: -1 south, 1 north.
    degacc : float, default=0.01
        Accuracy (default 0.01 deg) for numerical discretization of the Green's functions.
    requested_outputs : list, default=['default']
        Additional outputs requested (default: EsaUmotion).
    transitions : list, default=[]
        Indices into parts of the mesh that will be icecaps.

    Methods
    -------
    __init__(self, other=None)
        Initializes the esa parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the esa parameters.
    __str__(self)
        Returns a short string identifying the class.
    process_outputs(self, md=None, return_default_outputs=False)
        Process requested outputs, expanding 'default' to appropriate outputs.
    marshall_class(self, fid, prefix, md=None)
        Marshall parameters to a binary file.

    Examples
    --------
    md.esa = pyissm.param.esa()
    md.esa.deltathickness = thickness_change
    md.esa.love_h = 0.6
    md.esa.love_l = 0.1
    md.esa.degacc = 0.005
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

        s += '{}\n'.format(param_utils.fielddisplay(self, 'deltathickness', 'thickness change: ice height equivalent [m]'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'love_h', 'load Love number for radial displacement'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'love_l', 'load Love number for horizontal displaements'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'hemisphere', 'North-south, East-west components of 2-D horiz displacement vector:-1 south, 1 north'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'degacc', 'accuracy (default .01 deg) for numerical discretization of the Green''s functions'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'transitions', 'indices into parts of the mesh that will be icecaps'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'requested_outputs', 'additional outputs requested (default: EsaUmotion)'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - esa Class'
        return s

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
                
        ## Write fields
        execute.WriteData(fid, prefix, obj = self, fieldname = 'deltathickness', format = 'DoubleMat', mattype = 2)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'love_h', format = 'DoubleMat', mattype = 1)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'love_l', format = 'DoubleMat', mattype = 1)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'hemisphere', format = 'Integer')
        execute.WriteData(fid, prefix, obj = self, fieldname = 'degacc', format = 'Double')
        execute.WriteData(fid, prefix, obj = self, fieldname = 'transitions', format = 'MatArray')
        execute.WriteData(fid, prefix, name = 'md.esa.requested_outputs', data = self.process_outputs(md), format = 'StringArray')