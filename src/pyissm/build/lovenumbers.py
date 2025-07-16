import numpy as np
from . import build_utils
from . import class_registry

@class_registry.register_class
class lovenumbers(class_registry.manage_state):
    '''
    lovenumbers Class definition
    '''

    # Initialise with default parameters
    def __init__(self, other = None):
        # Loading love numbers
        self.h = np.nan # Provided by PREM model
        self.k = np.nan # idem
        self.l = np.nan # idem

        # Tidal love numbers for computing rotational feedback
        self.th = np.nan
        self.tk = np.nan
        self.tl = np.nan
        self.tk2secular = 0 # deg 2 secular number
        self.pmtf_colinear = np.nan
        self.pmtf_ortho = np.nan

        # Time/frequency for visco-elastic love numbers
        self.timefreq = np.nan
        self.istime = 1

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '---------------------------------------\n'
        s += '****      NOT YET IMPLEMENTED      ****\n'
        s += '---------------------------------------\n\n'
        s += '   lovenumbers parameters:\n'
        s += '{}\n'.format(build_utils.fielddisplay(self, 'h', 'load Love number for radial displacement'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'k', 'load Love number for gravitational potential perturbation'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'l', 'load Love number for horizontal displacements'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'th', 'tidal load Love number (deg 2)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'tk', 'tidal load Love number (deg 2)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'tl', 'tidal load Love number (deg 2)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'tk2secular', 'secular fluid Love number'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'pmtf_colinear', 'Colinear component of the Polar Motion Transfer Function (e.g. x-motion due to x-component perturbation of the inertia tensor)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'pmtf_ortho', 'Orthogonal component of the Polar Motion Transfer Function (couples x and y components, only used for Chandler Wobble)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'istime', 'time (default: 1) or frequency love numbers (0)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'timefreq', 'time/frequency vector (yr or 1/yr)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'pmtf_colinear', 'Colinear component of the Polar Motion Transfer Function (e.g. x-motion due to x-component perturbation of the inertia tensor)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'pmtf_ortho', 'Orthogonal component of the Polar Motion Transfer Function (couples x and y components, only used for Chandler Wobble)'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - lovenumbers Class'
        return s

