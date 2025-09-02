from . import param_utils
from . import class_registry
import numpy as np
from .. import execute

@class_registry.register_class
class rifts(class_registry.manage_state):
    """
    Rifts parameters class for ISSM.

    This class encapsulates parameters for modeling rifts in the ISSM (Ice Sheet System Model) framework.
    Rifts are fractures or cracks in ice sheets that can affect ice dynamics and calving processes.
    This class stores structural information about rifts including their geometry, properties, and melange characteristics.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in `other` differ from default values, they will override the default values.

    Attributes
    ----------
    riftstruct : list, default=[]
        Structure containing all rift information (vertices coordinates, segments, type of melange, etc.).
    riftproperties : list, default=[]
        Rift properties including physical and mechanical characteristics.

    Methods
    -------
    __init__(self, other=None)
        Initializes the rifts parameters, optionally inheriting from another instance.
    __repr__(self)
        Returns a detailed string representation of the rifts parameters.
    __str__(self)
        Returns a short string identifying the class.
    marshall_class(self, fid, prefix, md=None)
        Marshall parameters to a binary file

    Examples
    --------
    md.rifts = pyissm.param.rifts()
    md.rifts.riftstruct = rift_structure_data
    md.rifts.riftproperties = rift_properties_data
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.riftstruct = []
        self.riftproperties = []

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   rift parameters:\n'

        s += '{}\n'.format(param_utils.fielddisplay(self, 'riftstruct', 'structure containing all rift information (vertices coordinates, segments, type of melange, ...)'))
        s += '{}\n'.format(param_utils.fielddisplay(self, 'riftproperties', 'rift properties'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - rifts Class'
        return s


    # Marshall method for saving the rifts parameters
    def marshall_class(self, fid, prefix, md = None):
        """
        Marshall [rifts] parameters to a binary file.

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

        ## Process rift information
        if (not self.riftstruct) or (
            not isinstance(self.riftstruct, (tuple, list, dict))
            and np.any(np.isnan(self.riftstruct))
            ):
            numrifts = 0
        else:
            numrifts = len(self.riftstruct)

        numpairs = 0

        if numrifts > 0:
            for rift in self.riftstruct:
                numpairs += np.size(rift['penaltypairs'], axis = 0)

            ## Convert strings in riftstruct to hard coded numbers:
            FillDict = {'Air': 0,
                        'Ice': 1,
                        'Melange': 2,
                        'Water': 3}
            
            for rift in self.riftstruct:
                if rift['fill'] in ['Air', 'Ice', 'Melange', 'Water']:
                    rift['fill'] = FillDict[rift['fill']]

            # +2 for nodes + 2 for elements + 2 for  normals + 1 for length + 1 for fill + 1 for friction + 1 for fraction + 1 for fractionincrement + 1 for state.
            data = np.zeros((numpairs, 12))
            count = 0
            for rift in self.riftstruct:
                numpairsforthisrift = np.size(rift['penaltypairs'], 0)
                data[count:count + numpairsforthisrift, 0:7] = rift['penaltypairs']
                data[count:count + numpairsforthisrift, 7] = rift['fill']
                data[count:count + numpairsforthisrift, 8] = rift['friction']
                data[count:count + numpairsforthisrift, 9] = rift['fraction']
                data[count:count + numpairsforthisrift, 10] = rift['fractionincrement']
                data[count:count + numpairsforthisrift, 11] = rift['state'].reshape(-1)
                count += numpairsforthisrift
            else:
                data = np.zeros((numpairs, 12))
        
        ## Write fields
        execute.WriteData(fid, prefix, name = 'md.rifts.numrifts', data = numrifts, format = 'Integer')
        execute.WriteData(fid, prefix, name = 'md.rifts.riftstruct', data = data, format = 'DoubleMat', mattype = 3)