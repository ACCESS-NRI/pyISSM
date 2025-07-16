import numpy as np

CLASS_REGISTRY = {}
def register_class(cls):
    parts = cls.__module__.split(".")
    if "build" in parts:
        build_index = parts.index("build")
        key_parts = parts[build_index + 1:]  # skip 'build'
    else:
        key_parts = parts[-2:]  # fallback

    key = ".".join(key_parts + [cls.__name__])
    CLASS_REGISTRY[key] = cls
    return cls

def map_classtype(classtype: str) -> str:
    '''
    Masp legacy class types to modern ones for backwards compatability.
    '''
    legacy_class_map = {
        'basalforcings.basalforcings': 'basalforcings.default',
        'basalforcingspico.basalforcingspico': 'basalforcings.pico',
        'linearbasalforcings.linearbasalforcings': 'basalforcings.linear',
        'linearbasalforcingsarma.linearbasalforcingsarma': 'basalforcings.lineararma',
        'mismipbasalforcings.mismipbasalforcings': 'basalforcings.mismip',
        'plumebasalforcings.plumebasalforcings': 'basalforcings.plume',
        'calving.calving': 'calving.default',
        'calvingcrevassedepth.calvingcrevassedepth': 'calving.crevassedepth',
        'calvingdev.calvingdev': 'calving.dev',
        'calvinglevermann.calvinglevermann': 'calving.levermann',
        'calvingparameterization.calvingparameterization': 'calving.parameterization',
        'calvingvonmises.calvingvonmises': 'calving.vonmises',
        'dsl.dsl': 'dsl.default',
        'dslmme.dslmme': 'dsl.mme',
        'love.love': 'love.default',
        'fourierlove': 'love.fourier',
        'friction.friction': 'friction.default',
        'frictioncoulomb.frictioncoulomb': 'friction.coulomb',
        'frictioncoulomb2.frictioncoulomb2': 'friction.coulomb2',
        'frictionhydro.frictionhydro': 'friction.hydro',
        'frictionjosh.frictionjosh': 'friction.josh',
        'frictionpism.frictionpism': 'friction.pism',
        'frictionregcoulomb.frictionregcoulomb': 'friction.regcoulomb',
        'frictionregcoulomb2.frictionregcoulomb2': 'friction.regcoulomb2',
        'frictionschoof.frictionschoof': 'friction.schoof',
        'frictionshakti.frictionshakti': 'friction.shakti',
        'frictionwaterlayer.frictionwaterlayer': 'friction.waterlayer',
        'frictionweertman.frictionweertman': 'friction.weertman',
        'frontalforcings.frontalforcings': 'frontalforcings.default',
        'frontalforcingsrignot.frontalforcingsrignot': 'frontalforcings.rignot',
        'frontalforcingsrignotarma.frontalforcingsrignotarma': 'frontalforcings.rignotarma',
        'hydrologyarmapw.hydrologyarmapw': 'hydrology.armapw',
        'hydrologydc.hydrologydc': 'hydrology.dc',
        'hydrologyglads.hydrologyglads': 'hydrology.glads',
        'hydrologypism.hydrologypism': 'hydrology.pism',
        'hydrologyshakti.hydrologyshakti': 'hydrology.shakti',
        'hydrologyshreve.hydrologyshreve': 'hydrology.shreve',
        'hydrologytws.hydrologytws': 'hydrology.tws',
        'inversion.inversion': 'inversion.default',
        'taoinversion.taoinversion': 'inversion.tao',
        'm1qn3inversion.m1qn3inversion': 'inversion.m1qn3',
        'matice.matice': 'materials.ice',
        'matdamageice.matdamageice': 'materials.damageice',
        'matenhancedice.matenhancedice': 'materials.enhancedice',
        'matestar.matestar': 'materials.estar',
        'mesh2d.mesh2d': 'mesh.mesh2d',
        'mesh3d.mesh3d': 'mesh.mesh3d',
        'mesh3dprisms.mesh3dprisms': 'mesh.mesh3dprisms',
        'mesh2dvertical.mesh2dvertical': 'mesh.mesh2dvertical',
        'mesh3dsurface.mesh3dsurface': 'mesh.mesh3dsurface',
        'SMBforcing.SMBforcing': 'smb.default',
        'SMBarma.SMBarma': 'smb.arma',
        'SMBcomponents.SMBcomponents': 'smb.components',
        'SMBd18opdd.SMBd18opdd': 'smb.d18opdd',
        'SMBgradients.SMBgradients': 'smb.gradients',
        'SMBgradientscomponents.SMBgradientscomponents': 'smb.gradientscomponents',
        'SMBgradientsela.SMBgradientsela': 'smb.gradientsela',
        'SMBhenning.SMBhenning': 'smb.henning',
        'SMBmeltcomponents.SMBmeltcomponents': 'smb.meltcomponents',
        'SMBpdd.SMBpdd': 'smb.pdd',
        'SMBpddSicopolis.SMBpddSicopolis': 'smb.pddSicopolis',
        'solidearth.solidearth': 'solidearth.earth', # TODO: Check for earth or europa?
        'solidearthsettings.solidearthsettings': 'solidearth.settings',
        'solidearthsolution.solidearthsolution': 'solidearth.solution',
        'timestepping.timestepping': 'timestepping.default',
        'timesteppingadaptive.timesteppingadaptive': 'timestepping.adaptive'
    }
    if classtype in legacy_class_map:
        print(f"⚠️ Legacy classtype '{classtype}' mapped to '{legacy_class_map[classtype]}'")
    return legacy_class_map.get(classtype, classtype)


def create_instance(classtype: str):
    classtype = map_classtype(classtype)
    if classtype not in CLASS_REGISTRY:
        print(f"⚠️ Unknown classtype {classtype}. Skipping...")
        return None
        # raise ValueError(f'create_instance: Unknown class type: {classtype}')
    return CLASS_REGISTRY[classtype]()

## Manage state for save/load and inheritance
class manage_state:

    ## Allow inheritance from existing model class
    def __init__(self, other=None):

        # If other is provided...
        if other is not None:
            # Loop through all attributes of the current class...
            for attr in vars(self):
                # If the same attribute existing in the provided class, get the two fields
                if hasattr(other, attr):
                    field_other = getattr(other, attr)
                    field_self = getattr(self, attr)

                    # If the fields are different, replace the field_self with field_other
                    if not self._fields_equal(field_self, field_other):
                        setattr(self, attr, field_other)

    ## Check if fields are equal (used above)
    def _fields_equal(self, a, b):

        # If both values are NaN...
        if isinstance(a, float) and np.isnan(a):
            return isinstance(b, float) and np.isnan(b)
        if isinstance(b, float) and np.isnan(b):
            return isinstance(a, float) and np.isnan(a)

        # If both values are equal arrays...
        if isinstance(a, np.ndarray) and isinstance(b, np.ndarray):
            return np.array_equal(a, b)

        # If both values are equal scalars...
        return a == b

    ## Get the current state of self
    def __getstate__(self):
        return self.__dict__

    ## Set the current stat of self
    def __setstate__(self, state):
        self.__dict__.update(state)