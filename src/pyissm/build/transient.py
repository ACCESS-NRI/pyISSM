from . import build_utils
from . import class_registry

@class_registry.register_class
class transient(class_registry.manage_state):
    """
    Transient solution parameters class for ISSM.

    This class encapsulates parameters for configuring transient (time-dependent) simulations in the ISSM (Ice Sheet System Model) framework.
    It allows users to enable or disable various physics components and models that can be included in transient simulations,
    such as age tracking, surface mass balance, thermal evolution, and grounding line migration.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    isage : int, default=0
        Indicates if age model is requested in the transient.
    issmb : int, default=1
        Indicates if a surface mass balance solution is used in the transient.
    ismasstransport : int, default=1
        Indicates if a masstransport solution is used in the transient.
    ismmemasstransport : int, default=0
        Indicates whether an MME masstransport solution is used in the transient.
    isoceantransport : int, default=0
        Indicates whether an ocean masstransport solution is used in the transient.
    isstressbalance : int, default=1
        Indicates if a stressbalance solution is used in the transient.
    isthermal : int, default=1
        Indicates if a thermal solution is used in the transient.
    isgroundingline : int, default=0
        Indicates if a groundingline migration is used in the transient.
    isesa : int, default=0
        Indicates whether an elastic adjustment model is used in the transient.
    isdamageevolution : int, default=0
        Indicates whether damage evolution is used in the transient.
    ismovingfront : int, default=0
        Indicates whether a moving front capability is used in the transient.
    ishydrology : int, default=0
        Indicates whether an hydrology model is used.
    isdebris : int, default=0
        Indicates whether a debris model is used.
    issampling : int, default=0
        Indicates whether sampling is used in the transient.
    isslc : int, default=0
        Indicates if a sea level change solution is used in the transient.
    amr_frequency : int, default=0
        Frequency at which mesh is refined in simulations with multiple time_steps.
    isoceancoupling : int, default=0
        Indicates whether coupling with an ocean model is used in the transient (1 for cartesian coordinates, 2 for lat/long coordinates).
    requested_outputs : str, default='List of requested outputs'
        List of additional outputs requested.

    Methods
    -------
    __init__(self, other=None)
        Initializes the transient parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the transient parameters.
    __str__(self)
        Returns a short string identifying the class.

    Examples
    --------
    md.transient = pyissm.build.transient()
    md.transient.isage = 1
    md.transient.isgroundingline = 1
    md.transient.requested_outputs = ['IceVolume', 'IceVolumeAboveFloatation']
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.isage = 0
        self.issmb = 1
        self.ismasstransport = 1
        self.ismmemasstransport = 0
        self.isoceantransport = 0
        self.isstressbalance = 1
        self.isthermal = 1
        self.isgroundingline = 0
        self.isesa = 0
        self.isdamageevolution = 0
        self.ismovingfront = 0
        self.ishydrology = 0
        self.isdebris = 0
        self.issampling = 0
        self.isslc = 0
        self.amr_frequency = 0
        self.isoceancoupling = 0
        self.requested_outputs = 'List of requested outputs'

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   transient solution parameters:\n'

        s += '{}\n'.format(build_utils.fielddisplay(self, 'isage', 'indicates if age model is requested in the transient'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'issmb', 'indicates if a surface mass balance solution is used in the transient'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'ismasstransport', 'indicates if a masstransport solution is used in the transient'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'ismmemasstransport', 'indicates whether an MME masstransport solution is used in the transient'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'isoceantransport', 'indicates whether an ocean masstransport solution is used in the transient'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'isstressbalance', 'indicates if a stressbalance solution is used in the transient'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'isthermal', 'indicates if a thermal solution is used in the transient'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'isgroundingline', 'indicates if a groundingline migration is used in the transient'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'isesa', 'indicates whether an elastic adjustment model is used in the transient'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'isdamageevolution', 'indicates whether damage evolution is used in the transient'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'ismovingfront', 'indicates whether a moving front capability is used in the transient'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'ishydrology', 'indicates whether an hydrology model is used'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'isdebris', 'indicates whether a debris model is used'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'issampling', 'indicates whether sampling is used in the transient'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'isslc', 'indicates if a sea level change solution is used in the transient'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'isoceancoupling', 'indicates whether coupling with an ocean model is used in the transient (1 for cartesian coordinates, 2 for lat/long coordinates'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'amr_frequency', 'frequency at which mesh is refined in simulations with multiple time_steps'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'requested_outputs', 'list of additional outputs requested'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - transient Class'
        return s

