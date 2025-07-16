from. import build_utils
from . import class_registry

@class_registry.register_class
class issmsettings(class_registry.manage_state):
    '''
    issmsettings Class definition
    '''

    # Initialise with default parameters
    def __init__(self, other = None):
        self.results_on_nodes = 'List of output'
        self.io_gather = 1
        self.lowmem = 0
        self.output_frequency = 1
        self.sb_coupling_frequency = 1
        self.checkpoint_frequency = 0
        self.waitonlock = pow(2, 31) - 1
        self.solver_residue_threshold = 1e-6

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   general issmsettings parameters:\n'

        s += '{}\n'.format(build_utils.fielddisplay(self, 'results_on_nodes', "list of output for which results will be output for all the nodes of each element, Use 'all' for all output on nodes."))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'io_gather', 'I / O gathering strategy for result outputs (default 1)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'lowmem', 'is the memory limited ? (0 or 1)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'output_frequency', 'number of time steps between two saves (e.g., 5 means that results are only saved every 5 time steps)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'sb_coupling_frequency', 'frequency at which StressBalance solver is coupled (default 1)'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'checkpoint_frequency', 'frequency at which the runs are being recorded, allowing for a restart'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'waitonlock', 'maximum number of minutes to wait for batch results, or return 0'))
        s += '{}\n'.format(build_utils.fielddisplay(self, 'solver_residue_threshold', 'throw an error if solver residue exceeds this value (NaN to deactivate)'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - issmsettings Class'
        return s

