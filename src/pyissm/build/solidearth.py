from .. import utils
from . import build_utils
from . import class_registry
from . import rotational
from . import lovenumbers

## ------------------------------------------------------
## solidearth.earth
## ------------------------------------------------------
@class_registry.register_class
class earth(class_registry.manage_state):
    '''
    solidearth.earth Class definition
    '''

    # Initialise with default parameters
    def __init__(self, other = None):
        self.settings          = settings()
        self.external          = None
        self.lovenumbers       = lovenumbers()
        self.rotational        = rotational()
        self.planetradius      = utils.planetradius('earth')
        self.requested_outputs = 'List of requested outputs'
        self.transfercount     = 'List of transfer count'
        self.transitions       = 'List of transitions'
        self.partitionice      = 'List of partionice'
        self.partitionhydro    = 'List of partitionhydro'
        self.partitionocean    = 'List of partitionocean'

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   solidearth inputs, forcings, and settings:\n'
        s += '{}\n'.format(build_utils.fielddisplay(self, 'planetradius', 'planet radius [m]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'transitions', 'indices into parts of the mesh that will be icecaps'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'transfercount', 'number of icecaps vertices are part of'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'requested_outputs', 'additional outputs requested'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'partitionice', 'ice partition vector for barystatic contribution'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'partitionhydro', 'hydro partition vector for barystatic contribution'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'partitionocean', 'ocean partition vector for barystatic contribution'))
        if not self.external:
            s += '{}\n'.format(build_utils.fielddisplay(self, 'external', 'external solution, of the type solidearthsolution'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - solidearth.earth Class'
        return s

## ------------------------------------------------------
## solidearth.europa
## ------------------------------------------------------
@class_registry.register_class
class europa(class_registry.manage_state):
    '''
    solidearth.europa Class definition
    '''

    # Initialise with default parameters
    def __init__(self, other = None):
        self.settings          = settings()
        self.external          = None
        self.lovenumbers       = lovenumbers()
        self.rotational        = rotational()
        self.planetradius      = utils.planetradius('europa')
        self.requested_outputs = 'List of requested outputs'
        self.transfercount     = 'List of transfer count'
        self.transitions       = 'List of transitions'
        self.partitionice      = 'List of partionice'
        self.partitionhydro    = 'List of partitionhydro'
        self.partitionocean    = 'List of partitionocean'

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   solidearth inputs, forcings, and settings:\n'
        s += '{}\n'.format(build_utils.fielddisplay(self, 'planetradius', 'planet radius [m]'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'transitions', 'indices into parts of the mesh that will be icecaps'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'transfercount', 'number of icecaps vertices are part of'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'requested_outputs', 'additional outputs requested'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'partitionice', 'ice partition vector for barystatic contribution'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'partitionhydro', 'hydro partition vector for barystatic contribution'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'partitionocean', 'ocean partition vector for barystatic contribution'))
        if not self.external:
            s += '{}\n'.format(build_utils.fielddisplay(self, 'external', 'external solution, of the type solidearthsolution'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - solidearth.europa Class'
        return s

## ------------------------------------------------------
## solidearth.settings
## ------------------------------------------------------
@class_registry.register_class
class settings(class_registry.manage_state):
    '''
    solidearth.settings Class definition
    '''

    # Initialise with default parameters
    def __init__(self, other = None):
        self.reltol = 0
        self.abstol = 0
        self.maxiter = 0
        self.selfattraction = 1
        self.elastic = 1
        self.viscous = 1
        self.rotation = 1
        self.grdocean = 1
        self.ocean_area_scaling = 0
        self.runfrequency = 1
        self.sealevelloading = 1
        self.isgrd = 0
        self.compute_bp_grd = 0
        self.degacc = 0
        self.timeacc = 0
        self.horiz = 0
        self.grdmodel = 0
        self.cross_section_shape = 0

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   solidearth settings:\n'

        s += '{}\n'.format(build_utils.fielddisplay(self, 'reltol', 'sea level change relative convergence criterion (default, NaN: not applied)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'abstol', 'sea level change absolute convergence criterion (default, NaN: not applied)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'maxiter', 'maximum number of nonlinear iterations'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'grdocean', 'does this planet have an ocean, if set to 1: global water mass is conserved in GRD module (default: 1)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'ocean_area_scaling', 'correction for model representation of ocean area (default: No correction)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'sealevelloading', 'enables surface loading from sea-level change (default: 1)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'isgrd', 'compute GRD patterns (default: 1'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'compute_bp_grd', 'compute GRD patterns for bottom pressure loads (default 1)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'runfrequency', 'how many time steps we skip before we run solidearthsettings solver during transient (default: 1)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'selfattraction', 'enables surface mass load to perturb the gravity field'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'elastic', 'enables elastic deformation from surface loading'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'viscous', 'enables viscous deformation from surface loading'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'rotation', 'enables polar motion to feedback on the GRD fields'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'degacc', 'accuracy (default: .01 deg) for numerical discretization of the Green\'s functions'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'timeacc', 'time accuracy (default: 1 year) for numerical discretization of the Green\'s functions'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'grdmodel', 'type of deformation model, 0 for no GRD, 1 for spherical GRD model (SESAW model), 2 for half-space planar GRD (visco-elastic model from Ivins)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'cross_section_shape', '1: square-edged (default). 2: elliptical. See iedge in GiaDeflectionCore'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - solidearth.settings Class'
        return s

## ------------------------------------------------------
## solidearth.solution
## ------------------------------------------------------
@class_registry.register_class
class solution(class_registry.manage_state):
    '''
    solidearth.solution Class definition
    '''

    # Initialise with default parameters
    def __init__(self, other = None):
        self.displacementeast = 'displacementeast timeseries'
        self.displacementnorth = 'displacementsnorth timeseries'
        self.displacementup = 'displacementup timeseries'
        self.geoid = 'geoid timeseries'

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   solidearth solution:\n'
        s += '         units for time series is (yr)\n'

        s += '{}\n'.format(build_utils.fielddisplay(self, 'displacementeast', 'solid-Earth Eastwards bedrock displacement series (m)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'displacementnorth', 'solid-Earth Northwards bedrock displacement time series (m)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'displacementup', 'solid-Earth bedrock uplift time series (m)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'geoid', 'solid-Earth geoid time series (m)'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - solidearth.solution Class'
        return s