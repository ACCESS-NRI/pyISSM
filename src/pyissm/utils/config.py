"""
Configuration functions for ISSM

This module contains various functions that relate to the configuration of the ISSM codebase.
"""

import collections
from . import wrappers

def mumps_options(**kwargs):
    """
    Set the MUMPS options
    """

    petsc_major = wrappers.IssmConfig('_PETSC_MAJOR_')[0]
    petsc_minor = wrappers.IssmConfig('_PETSC_MINOR_')[0]

    ## Define defaults (conditional on petsc version)
    if petsc_major == 2:
        defaults = {'toolkit': 'petsc',
                    'mat_type': 'aijumps',
                    'ksp_type': 'preonly',
                    'pc_type': 'lu',
                    'mat_mumps_icntl_14': 120}
    
    if petsc_major == 3:
        defaults = {'toolkit': 'petsc',
            'mat_type': 'mpiaij',
            'ksp_type': 'preonly',
            'pc_type': 'lu',
            'mat_mumps_icntl_14': 120}
        if petsc_minor > 8:
            defaults['pc_factor_mat_solver_type'] = 'mumps'
        else:
            defaults['pc_factor_mat_solver_package'] = 'mumps'
        
        defaults['mat_mumps_icntl_28'] = 1 # 1 = serial; 2 = parralel
        defaults['mat_mumps_icntl_29'] = 2 # parallel ordering: 1 = ptscotch, 2 = parmetis

    ## Update with user options
    mumps = collections.OrderedDict(defaults)
    mumps.update(kwargs)
    return mumps

def iluasm_options(**kwargs):
    """
    ILUASMOPTIONS -

       Usage:
          options = iluasmoptions
    """

    ## Define defaults
    defaults = {'mat_type': 'aij',
                'ksp_type': 'gmres',
                'pc_type': 'asm',
                'sub_pc_type': 'ilu',
                'pc_asm_overlap': 5,
                'ksp_max_it': 100,
                'ksp_rtol': 1e-15}

    ## Update with user options
    iluasm = collections.OrderedDict(defaults)
    iluasm.update(kwargs)

    return iluasm

def issm_mumps_solver(**kwargs):
    #ISSMSOLVE - return issm solver options
    #
    #   Usage:
    #      options = issmsolver

    ## Define defaults
    defaults = {'toolkit': 'issm',
                'mat_type': 'mpiparse',
                'vec_type': 'mpi',
                'solver_type': 'mumps'}
    
    ## Update with user options
    mumps = collections.OrderedDict(defaults)
    mumps.update(kwargs)
    
    return mumps

def issm_gsl_solver(**kwargs):
    #ISSMSOLVE - return issm solver options
    #
    #   Usage:
    #      options = issmsolver

    ## Define defaults
    defaults = {'toolkit': 'issm',
                'mat_type': 'dense',
                'vec_type': 'seq',
                'solver_type': 'gsl'}
    
    ## Update with user options
    gsl = collections.OrderedDict(defaults)
    gsl.update(kwargs)

    return gsl
