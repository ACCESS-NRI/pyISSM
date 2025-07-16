from . import build_utils
from . import class_registry

## ------------------------------------------------------
## results.default
## ------------------------------------------------------
@class_registry.register_class
class default(class_registry.manage_state):

    # Initialise with default parameters
    def __init__(self, other = None):
        pass

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):  #{{{
        s = ''
        for key, value in self.__dict__.items():
            if isinstance(value, resultsdakota):
                lengthvalue = 1
            else:
                try:
                    lengthvalue = len(value)
                except TypeError:
                    lengthvalue = 1
            s += '    {}: [1x{} struct]\n'.format(key, lengthvalue)
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - results.default Class'
        return s

## ------------------------------------------------------
## results.resultsdakota
## ------------------------------------------------------
@class_registry.register_class
class resultsdakota(class_registry.manage_state):

    # Initialise with default parameters
    def __init__(self, other = None):
        pass

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = ''
        for key, value in self.__dict__.items():
            s += '    {}: '.format(key)
            if isinstance(value, list):
                s += '[{} element list]'.format(len(value))
            else:
                s += '{}'.format(value)
            s += '\n'
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - resultsdakota Class'
        return s

## ------------------------------------------------------
## results.solution
## ------------------------------------------------------
@class_registry.register_class
class solution(class_registry.manage_state):

    # Initialise with default parameters
    def __init__(self, *args):
        self.steps = None
        if len(args) == 1:
            arg = args[0]
            if isinstance(arg, list):
                self.steps = arg
            else:
                raise Exception('solution class error: if initializing with an argument, that argument should be an empty list or a list of instances of solutionstep')
        else:
            self.steps = [solutionstep()]

    # Define repr
    def __repr__(self):  #{{{
        s = ''
        numsteps = len(self.steps)
        if numsteps == 1:
            for key, value in self.steps[0].__dict__.items():
                s += '    {}: {}\n'.format(key, value)
        else:
            s = '  1x{} struct array with fields:\n'.format(numsteps)
            s += '\n'
            for fieldname in self.steps[0].getfieldnames():
                s += '    {}\n'.format(fieldname)
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - solution Class'
        return s

## ------------------------------------------------------
## results.solutionstep
## ------------------------------------------------------
@class_registry.register_class
class solutionstep(class_registry.manage_state):

    # Initialise with default parameters
    def __init__(self, other = None):
        pass

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = ''
        width = build_utils.getlongestfieldname(self)
        for key, value in self.__dict__.items():
            s += '    {:{width}s}: {}\n'.format(key, value, width=width)
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - solutionstep Class'
        return s