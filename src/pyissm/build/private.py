from . import build_utils
from . import class_registry

@class_registry.register_class
class private(class_registry.manage_state):
    '''
    private Class definition
    '''

    # Initialise with default parameters
    def __init__(self, other = None):
        self.isconsistent = True
        self.runtimename = ''
        self.bamg = 'OrderedDict()'
        self.solution = ''

        # Inherit matching fields from provided class
        super().__init__(other)

        ## TODO: How should OrderedDict() be implemented?

    # Define repr
    def __repr__(self):
        s = '   private parameters -- do not change:\n'

        s += '{}\n'.format(build_utils.fielddisplay(self, 'isconsistent', 'is model self consistent?'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'runtimename', 'name of the run launched'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'bamg', 'structure with mesh properties constructed if bamg is used to mesh the domain'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'solution', 'type of solution launched'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - private Class'
        return s

