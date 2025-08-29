import numpy as np
import collections
import warnings
from . import param_utils
from . import class_registry
from .. import execute

## ------------------------------------------------------
## qmu.default
## ------------------------------------------------------
@class_registry.register_class
class default(class_registry.manage_state):
    '''
    qmu.default Class definition
    '''

    # Initialise with default parameters
    def __init__(self, other = None):
        self.isdakota = 0
        self.output = 0
        self.variables = 'OrderedStruct() -- NOT IMPLEMENTED'
        self.correlation_matrix = 'List of correlation matrix'
        self.responses = 'OrderedStruct() -- NOT IMPLEMENTED'
        self.method = collections.OrderedDict()
        self.params = 'OrderedStruct() -- NOT IMPLEMENTED'
        self.statistics = statistics()
        self.results = collections.OrderedDict()
        self.numberofresponses = 0
        self.variabledescriptors = 'List of variable descriptors'
        self.variablepartitions = 'List of variable partitions'
        self.variablepartitions_npart = 'List of variable partitions (npart)'
        self.variablepartitions_nt = 'List of variable partitions (nt)'
        self.responsedescriptors = 'List of response descriptors'
        self.responsepartitions = 'List of response partitions'
        self.responsepartitions_npart = 'List of response partitions (npart)'
        self.responsepartitions_nt = 'List of response partitions (nt)'
        self.mass_flux_profile_directory = np.nan
        self.mass_flux_profiles = np.nan
        self.mass_flux_segments = [] # List of mass flux segments
        self.adjacency = np.nan
        self.vertex_weight = np.nan

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '---------------------------------------\n'
        s += '****      NOT YET IMPLEMENTED      ****\n'
        s += '---------------------------------------\n\n'
        s += '   qmu parameters:\n'
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - qmu.default Class'
        return s

    # Marshall method for saving the qmu parameters
    def marshall_class(self, fid, prefix, md = None):
        """
        Marshall [qmu] parameters to a binary file.

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

        ## Write fields (Turn off)
        execute.WriteData(fid, prefix, name = 'md.qmu.isdakota', data = False, format = 'Boolean')
        execute.WriteData(fid, prefix, name = 'md.qmu.mass_flux_segments_present', data = False, format = 'Boolean')

        warnings.warn('pyissm.param.qmu::qmu not yet implemented. Turning off qmu.')


## ------------------------------------------------------
## qmu.statistics
## ------------------------------------------------------
@class_registry.register_class
class statistics(class_registry.manage_state):
    '''
    qmu.statistics Class definition
    '''

    # Initialise with default parameters
    def __init__(self, other = None):
        self.nfiles_per_directory = 5
        self.ndirectories = 50
        self.method = [{}]
        self.method[0]['name'] = 'None'
        self.method[0]['fields'] = []
        self.method[0]['steps'] = []
        self.method[0]['nbins'] = np.nan
        self.method[0]['indices'] = []

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '---------------------------------------\n'
        s += '****      NOT YET IMPLEMENTED      ****\n'
        s += '---------------------------------------\n\n'
        s += 'qmustatistics: post-Dakota run processing of QMU statistics:\n'
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - qmu.statistics Class'
        return s

