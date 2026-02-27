import numpy as np
import warnings
from pyissm.model.classes import class_utils, class_registry, hydrology
from pyissm.model import execute, mesh

@class_registry.register_class
class initialization(class_registry.manage_state):
    """
    Initialization class for ISSM.

    This class contains initial field values for various physical quantities in the ISSM framework.
    It provides storage for initial conditions of velocity, pressure, temperature, and other state variables
    that are used to initialize ice sheet model simulations.

    Parameters
    ----------
    other : any, optional
        Any other class object that contains common fields to inherit from. If values in ``other`` differ from default
        values, they will override the default values.

    Attributes
    ----------
    vx : :class:`numpy.ndarray`, default=nan
        x component of velocity [m/yr].
    vy : :class:`numpy.ndarray`, default=nan
        y component of velocity [m/yr].
    vz : :class:`numpy.ndarray`, default=nan
        z component of velocity [m/yr].
    vel : :class:`numpy.ndarray`, default=nan
        velocity norm [m/yr].
    pressure : :class:`numpy.ndarray`, default=nan
        pressure [Pa].
    temperature : :class:`numpy.ndarray`, default=nan
        temperature [K].
    enthalpy : :class:`numpy.ndarray`, default=nan
        enthalpy [J].
    waterfraction : :class:`numpy.ndarray`, default=nan
        fraction of water in the ice.
    sediment_head : :class:`numpy.ndarray`, default=nan
        sediment water head of subglacial system [m].
    epl_head : :class:`numpy.ndarray`, default=nan
        epl water head of subglacial system [m].
    epl_thickness : :class:`numpy.ndarray`, default=nan
        thickness of the epl [m].
    watercolumn : :class:`numpy.ndarray`, default=nan
        thickness of subglacial water [m].
    hydraulic_potential : :class:`numpy.ndarray`, default=nan
        Hydraulic potential (for GlaDS) [Pa].
    channelarea : :class:`numpy.ndarray`, default=nan
        subglacial water channel area (for GlaDS) [m2].
    sealevel : :class:`numpy.ndarray`, default=nan
        sea level [m].
    bottompressure : :class:`numpy.ndarray`, default=nan
        bottom pressure [Pa].
    dsl : :class:`numpy.ndarray`, default=nan
        dynamic sea level [m].
    str : :class:`numpy.ndarray`, default=nan
        surface temperature rate [K/yr].
    sample : :class:`numpy.ndarray`, default=nan
        Realization of a Gaussian random field.
    debris : :class:`numpy.ndarray`, default=nan
        Surface debris layer [m].
    age : :class:`numpy.ndarray`, default=nan
        Initial age [yr].

    Examples
    --------
    .. code-block:: python
    
        >>> md.initialization = pyissm.model.classes.initialization()
        >>> md.initialization.vx = vx_initial
        >>> md.initialization.vy = vy_initial
        >>> md.initialization.temperature = temp_initial
    """

    # Initialise with default parameters
    def __init__(self, other = None):
        self.vx = np.nan
        self.vy = np.nan
        self.vz = np.nan
        self.vel = np.nan # NOTE: This is not used/written to file.
        self.pressure = np.nan
        self.temperature = np.nan
        self.enthalpy = np.nan
        self.waterfraction = np.nan
        self.sediment_head = np.nan
        self.epl_head = np.nan
        self.epl_thickness = np.nan
        self.watercolumn = np.nan
        self.hydraulic_potential = np.nan
        self.channelarea = np.nan
        self.sealevel = np.nan
        self.bottompressure = np.nan
        self.dsl = np.nan
        self.str = np.nan
        self.sample = np.nan
        self.debris = np.nan
        self.age = np.nan

        # Inherit matching fields from provided class
        super().__init__(other)

    # Define repr
    def __repr__(self):
        s = '   initial field values:\n'

        s += '{}\n'.format(class_utils._field_display(self, 'vx', 'x component of velocity [m/yr]'))
        s += '{}\n'.format(class_utils._field_display(self, 'vy', 'y component of velocity [m/yr]'))
        s += '{}\n'.format(class_utils._field_display(self, 'vz', 'z component of velocity [m/yr]'))
        s += '{}\n'.format(class_utils._field_display(self, 'vel', 'velocity norm [m/yr]'))
        s += '{}\n'.format(class_utils._field_display(self, 'pressure', 'pressure [Pa]'))
        s += '{}\n'.format(class_utils._field_display(self, 'temperature', 'temperature [K]'))
        s += '{}\n'.format(class_utils._field_display(self, 'enthalpy', 'enthalpy [J]'))
        s += '{}\n'.format(class_utils._field_display(self, 'waterfraction', 'fraction of water in the ice'))
        s += '{}\n'.format(class_utils._field_display(self, 'watercolumn', 'thickness of subglacial water [m]'))
        s += '{}\n'.format(class_utils._field_display(self, 'sediment_head', 'sediment water head of subglacial system [m]'))
        s += '{}\n'.format(class_utils._field_display(self, 'epl_head', 'epl water head of subglacial system [m]'))
        s += '{}\n'.format(class_utils._field_display(self, 'epl_thickness', 'thickness of the epl [m]'))
        s += '{}\n'.format(class_utils._field_display(self, 'hydraulic_potential', 'Hydraulic potential (for GlaDS) [Pa]'))
        s += '{}\n'.format(class_utils._field_display(self, 'channelarea', 'subglaciale water channel area (for GlaDS) [m2]'))
        s += '{}\n'.format(class_utils._field_display(self, 'sample', 'Realization of a Gaussian random field'))
        s += '{}\n'.format(class_utils._field_display(self, 'debris', 'Surface debris layer [m]'))
        s += '{}\n'.format(class_utils._field_display(self, 'age', 'Initial age [yr]'))
        return s

    # Define class string
    def __str__(self):
        s = 'ISSM - initialization Class'
        return s
    
    # Extrude to 3D mesh
    def _extrude(self, md):
        """
        Extrude initialization fields to 3D
        """
        self.vx = mesh._project_3d(md, vector = self.vx, type = 'node')
        self.vy = mesh._project_3d(md, vector = self.vy, type = 'node')
        self.vz = mesh._project_3d(md, vector = self.vz, type = 'node')
        self.vel = mesh._project_3d(md, vector = self.vel, type = 'node')
        self.temperature = mesh._project_3d(md, vector = self.temperature, type = 'node')
        self.enthalpy = mesh._project_3d(md, vector = self.enthalpy, type = 'node')
        self.waterfraction = mesh._project_3d(md, vector = self.waterfraction, type = 'node')
        self.watercolumn = mesh._project_3d(md, vector = self.watercolumn, type = 'node')
        self.sediment_head = mesh._project_3d(md, vector = self.sediment_head, type = 'node', layer = 1)
        self.epl_head = mesh._project_3d(md, vector = self.epl_head, type = 'node', layer = 1)
        self.epl_thickness = mesh._project_3d(md, vector = self.epl_thickness, type = 'node', layer = 1)
        self.sealevel = mesh._project_3d(md, vector = self.sealevel, type = 'node', layer = 1)
        self.bottompressure = mesh._project_3d(md, vector = self.bottompressure, type = 'node', layer = 1)
        self.dsl = mesh._project_3d(md, vector = self.dsl, type = 'node', layer = 1)
        self.str = mesh._project_3d(md, vector = self.str, type = 'node', layer = 1)
        self.debris = mesh._project_3d(md, vector = self.debris, type = 'node', layer = 1)
        self.age = mesh._project_3d(md, vector = self.age, type = 'node', layer = 1)

        # Lithostatic pressure by default
        if np.ndim(md.geometry.surface) == 2:
            warnings.warn('pyissm.model.classes.initialization._extrude: Reshaping md.geometry.surface for your convenience but you should fix it in your model set up')
            self.pressure = md.constants.g * md.materials.rho_ice * (md.geometry.surface.reshape(-1, 1) - md.mesh.z)
        else:
            self.pressure = md.constants.g * md.materials.rho_ice * (md.geometry.surface - md.mesh.z)
            
        return self

    # Check model consistency
    def check_consistency(self, md, solution, analyses):
        """
        Check consistency of the [initialization] parameters.

        Parameters
        ----------
        md : :class:`pyissm.model.Model`
            The model object to check.
        solution : :class:`pyissm.model.solution`
            The solution object to check.
        analyses : list of :class:`str`
            List of analyses to check consistency for.

        Returns 
        -------
        md : :class:`pyissm.model.Model`
            The model object with any consistency errors noted.
        """

        ## StressbalanceAnalysis
        if 'StressbalanceAnalysis' in analyses and not solution == 'TransientSolution' and not md.transient.isstressbalance:
            if not np.any(np.isnan(md.initialization.vx)) and not np.any(np.isnan(md.initialization.vy)):
                if np.size(md.initialization.vx) > 1 or np.size(md.initialization.vy) > 1:
                    class_utils.check_field(md, fieldname = "initialization.vx", allow_nan = False, allow_inf = False, size = (md.mesh.numberofvertices,))
                    class_utils.check_field(md, fieldname = "initialization.vy", allow_nan = False, allow_inf = False, size = (md.mesh.numberofvertices,))

        ## MasstransportAnalysis
        if 'MasstransportAnalysis' in analyses and not solution == 'TransientSolution' and not md.transient.ismasstransport:
            class_utils.check_field(md, fieldname = "initialization.vx", allow_nan = False, allow_inf = False, size = (md.mesh.numberofvertices,))
            class_utils.check_field(md, fieldname = "initialization.vy", allow_nan = False, allow_inf = False, size = (md.mesh.numberofvertices,))

        ## OceantransportAnalysis
        if 'OceantransportAnalysis' in analyses:
            if solution == 'TransientSolution' and md.transient.isslc and md.transient.isoceantransport:
                class_utils.check_field(md, fieldname = "initialization.bottompressure", allow_nan = False, allow_inf = False, size = (md.mesh.numberofvertices,))
                class_utils.check_field(md, fieldname = "initialization.dsl", allow_nan = False, allow_inf = False, size = (md.mesh.numberofvertices,))
                class_utils.check_field(md, fieldname = 'initialization.str', sclara = True, allow_nan = False, allow_inf = False)
        
        ## BalancethicknessAnalysis
        if 'BalancethicknessAnalysis' in analyses and solution == 'BalancethicknessSolution':
            class_utils.check_field(md, fieldname = 'initialization.vx', size = (md.mesh.numberofvertices,), allow_nan = False, allow_inf = False)
            class_utils.check_field(md, fieldname = 'initialization.vy', size = (md.mesh.numberofvertices,), allow_nan = False, allow_inf = False)
            # Triangle with zero velocity
            vx_sum = np.sum(np.abs(md.initialization.vx[md.mesh.elements - 1]), axis=1)
            vy_sum = np.sum(np.abs(md.initialization.vy[md.mesh.elements - 1]), axis=1)
            ice_mask = np.min(md.mask.ice_levelset[md.mesh.elements - 1], axis=1) < 0
            if np.any((vx_sum == 0) & (vy_sum == 0) & ice_mask):
                md.check_message('at least one triangle has all its vertices with a zero velocity')
        
        ## ThermalAnalysis
        if 'ThermalAnalysis' in analyses and not solution == 'TransientSolution' and not md.transient.isthermal:
            class_utils.check_field(md, fieldname = 'initialization.vx', allow_nan = False, allow_inf = False, size = (md.mesh.numberofvertices,))
            class_utils.check_field(md, fieldname = 'initialization.vy', allow_nan = False, allow_inf = False, size = (md.mesh.numberofvertices,))
            class_utils.check_field(md, fieldname = 'initialization.temperature', allow_nan = False, allow_inf = False, size = 'universal')
            if md.mesh.dimension() == 3:
                class_utils.check_field(md, fieldname = 'initialization.vz', size = (md.mesh.numberofvertices, ), allow_nan = False, allow_inf = False)
            class_utils.check_field(md, fieldname = 'initialization.pressure', size = (md.mesh.numberofvertices, ), allow_nan = False, allow_inf = False)
        
        ## EnthalpyAnalysis
        if 'EnthalpyAnalysis' in analyses and md.thermal.isenthalpy:
            class_utils.check_field(md, fieldname = 'initialization.waterfraction', size = (md.mesh.numberofvertices, ), ge = 0)
            class_utils.check_field(md, fieldname = 'initialization.watercolumn', size = (md.mesh.numberofvertices, ), ge = 0)
            # pos = np.nonzero(md.initialization.waterfraction > 0.)[0]
            # if(pos.size):
            #     class_utils.check_field(md, fieldname = 'delta Tpmp', field = np.absolute(md.initialization.temperature[pos] - (md.materials.meltingpoint - md.materials.beta * md.initialization.pressure[pos])), lt = 1e-11, message = 'set temperature to pressure melting point at locations with waterfraction > 0')
        
        ## HydrologyShreveAnalysis
        if 'HydrologyShreveAnalysis' in analyses:
            if isinstance(md.hydrology, hydrology.shreve):
                if (solution == 'TransientSolution' and md.transient.ishydrology) or solution == 'HydrologySolution':
                    class_utils.check_field(md, fieldname = 'initialization.watercolumn', size = (md.mesh.numberofvertices, ), allow_nan = False, allow_inf = False)
        
        ## HydrologyTwsAnalysis
        if 'HydrologyTwsAnalysis' in analyses:
            if isinstance(md.hydrology, hydrology.tws):
                class_utils.check_field(md, fieldname = 'initialization.watercolumn', size = (md.mesh.numberofvertices, ), allow_nan = False, allow_inf = False)
        
        ## SealevelchangeAnalysis
        if 'SealevelchangeAnalysis' in analyses:
            if solution == 'TransientSolution' and md.transient.isslc:
                class_utils.check_field(md, fieldname = 'initialization.sealevel', size = (md.mesh.numberofvertices, ), allow_nan = False, allow_inf = False)
        
        ## HydrologyGlaDSAnalysis
        if 'HydrologyGlaDSAnalysis' in analyses:
            if isinstance(md.hydrology, hydrology.glads):
                class_utils.check_field(md, fieldname = 'initialization.watercolumn', size = (md.mesh.numberofvertices, ), allow_nan = False, allow_inf = False)
                class_utils.check_field(md, fieldname = 'initialization.hydraulic_potential', size = (md.mesh.numberofvertices, ), allow_nan = False, allow_inf = False)
                class_utils.check_field(md, fieldname = 'initialization.channelarea', ge = 0, size = (md.mesh.numberofelements, ), allow_nan = False, allow_inf = False)
        
        ## HydrologyDC (Inefficient/Efficient)
        if 'HydrologyDCInefficientAnalysis' in analyses:
            if isinstance(md.hydrology, hydrology.dc):
                class_utils.check_field(md, fieldname = 'initialization.sediment_head', size = (md.mesh.numberofvertices, ), allow_nan = False, allow_inf = False)
        
        
        if 'HydrologyDCEfficientAnalysis' in analyses:
            if isinstance(md.hydrology, hydrology.dc):
                if md.hydrology.isefficientlayer:
                    class_utils.check_field(md, fieldname = 'initialization.epl_head', size = (md.mesh.numberofvertices, ), allow_nan = False, allow_inf = False)
                    class_utils.check_field(md, fieldname = 'initialization.epl_thickness', size = (md.mesh.numberofvertices, ), allow_nan = False, allow_inf = False)
        
        ## SamplingAnalysis
        if 'SamplingAnalysis' in analyses and not solution == 'TransientSolution' and not md.transient.issampling:
            if np.any(np.isnan(md.initialization.sample)):
                class_utils.check_field(md, fieldname = 'initialization.sample', size = (md.mesh.numberofvertices, ), allow_nan = False, allow_inf = False)
        
        ## DebrisAnalysis
        if 'DebrisAnalysis' in analyses:
            if not np.isnan(md.initialization.debris):
                if (solution == 'TransientSolution' and md.transient.ishydrology) or solution == 'HydrologySolution':
                    class_utils.check_field(md, fieldname = 'initialization.debris', size = (md.mesh.numberofvertices, ), allow_nan = False, allow_inf = False)
        
        ## AgeAnalysis
        if 'AgeAnalysis' in analyses:
            if not np.isnan(md.initialization.age):
                if (solution == 'TransientSolution' and md.transient.ishydrology) or solution == 'HydrologySolution':
                    class_utils.check_field(md, fieldname = 'initialization.age', size = (md.mesh.numberofvertices, ), allow_nan = False, allow_inf = False)
        return md

    # Marshall method for saving the initialization parameters
    def marshall_class(self, fid, prefix, md = None):
        """
        Marshall [initialization] parameters to a binary file.

        Parameters
        ----------
        fid : :class:`file object`
            The file object to write the binary data to.
        prefix : :class:`str`
            Prefix string used for data identification in the binary file.
        md : :class:`pyissm.model.Model`, optional
            ISSM model object needed in some cases.
            
        Returns
        -------
        None
        """
        
        ## Write DoubleMat fields (consistent mattype)
        fieldnames = ['pressure', 'bottompressure', 'str', 'dsl', 'temperature',
            'waterfraction', 'sediment_head', 'epl_head', 'epl_thickness',
            'watercolumn', 'channelarea', 'hydraulic_potential', 'sample', 'debris']
        for field in fieldnames:
            execute.WriteData(fid, prefix, obj = self, fieldname = field, format = 'DoubleMat', mattype = 1)

        # Write other fields
        execute.WriteData(fid, prefix, obj = self, fieldname = 'vx', format = 'DoubleMat', mattype = 1, scale = 1 / md.constants.yts, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'vy', format = 'DoubleMat', mattype = 1, scale = 1 / md.constants.yts, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'vz', format = 'DoubleMat', mattype = 1, scale = 1 / md.constants.yts)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'sealevel', format = 'DoubleMat', mattype = 1, timeserieslength = md.mesh.numberofvertices + 1, yts = md.constants.yts)
        execute.WriteData(fid, prefix, obj = self, fieldname = 'age', format = 'DoubleMat', mattype = 1, scale = md.constants.yts)

        # Write conditional fields
        if md.thermal.isenthalpy:
            if (np.size(self.enthalpy) <= 1):
                # Reconstruct enthalpy
                self.enthalpy = np.zeros(md.mesh.numberofvertices)
                tpmp = md.materials.meltingpoint - md.materials.beta * md.initialization.pressure
                pos = md.initialization.waterfraction > 0
                self.enthalpy[:] = md.materials.heatcapacity * (md.initialization.temperature - md.constants.referencetemperature)
                self.enthalpy[pos] = md.materials.heatcapacity * (tpmp[pos] - md.constants.referencetemperature) + md.materials.latentheat * md.initialization.waterfraction[pos]

            execute.WriteData(fid, prefix, name = 'md.initialization.enthalpy', data = self.enthalpy, format = 'DoubleMat', mattype = 1)