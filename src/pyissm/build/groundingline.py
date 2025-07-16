from . import build_utils
from . import class_registry

@class_registry.register_class
class groundingline(class_registry.manage_state):
    '''
    groundingline Class definition
    '''

    # Initialise with default parameters
    def __init__(self, other = None):
        self.migration = 'SubelementMigration'
        self.friction_interpolation = 'SubelementFriction1'
        self.melt_interpolation = 'NoMeltOnPartiallyFloating'
        self.intrusion_distance = 0
        self.requested_outputs = 'List of requested outputs'

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   grounding line migration parameters:\n'

        s += '{}\n'.format(build_utils.fielddisplay(self, 'migration', 'type of grounding line migration: \'SoftMigration\', \'SubelementMigration\', \'AggressiveMigration\', \'Contact\', \'None\''))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'migration', 'type of friction interpolation on partially floating elements: ''SubelementFriction1'', ''SubelementFriction2'', ''NoFrictionOnPartiallyFloating'''))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'migration', 'type of melt interpolation on partially floating elements: \'SubelementMelt1\', \'SubelementMelt2\', \'IntrusionMelt\', \'NoMeltOnPartiallyFloating\', \'FullMeltOnPartiallyFloating\''))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'requested_outputs', 'additional outputs requested'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - groundingline Class'
        return s

