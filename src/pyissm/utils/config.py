"""
Configuration functions for ISSM

This module contains various functions that relate to the configuration of the ISSM codebase.
"""

import collections
import platform
import os
import socket
from pyissm.utils import wrappers

def mumps_options(**kwargs):
    """
    Set the MUMPS options for PETSc linear solver configuration.

    This function configures MUMPS (MUltifrontal Massively Parallel sparse direct Solver)
    options based on the PETSc version detected in the system. It provides version-specific
    defaults and allows user customization through keyword arguments.

    Parameters
    ----------
    **kwargs : dict
        Additional MUMPS options to override or supplement the defaults.

    Returns
    -------
    collections.OrderedDict
        Dictionary containing MUMPS configuration options with PETSc version-specific
        defaults updated with any user-provided options.

    Notes
    -----
    The function automatically detects PETSc major and minor versions and sets
    appropriate defaults:
    - PETSc v2.x: Uses 'aijumps' matrix type
    - PETSc v3.x: Uses 'mpiaij' matrix type with additional solver package configuration
    For PETSc v3.9+, uses 'pc_factor_mat_solver_type', otherwise uses the deprecated
    'pc_factor_mat_solver_package'.

    Examples
    --------
    >>> opts = mumps_options()
    >>> opts = mumps_options(mat_mumps_icntl_14=200)
    >>> opts = mumps_options(mat_mumps_icntl_28=2, mat_mumps_icntl_29=1)
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
    Set the ILU-ASM options for PETSc linear solver configuration.

    This function configures ILU (Incomplete LU factorization) with ASM (Additive Schwarz Method)
    preconditioner options for PETSc. It provides default settings for iterative solving
    and allows user customization through keyword arguments.

    Parameters
    ----------
    **kwargs : dict
        Additional ILU-ASM options to override or supplement the defaults.

    Returns
    -------
    collections.OrderedDict
        Dictionary containing ILU-ASM configuration options with defaults
        updated with any user-provided options.

    Notes
    -----
    The function sets up an iterative solver configuration with:
    - AIJ matrix type for sparse matrices
    - GMRES Krylov subspace method
    - Additive Schwarz Method (ASM) preconditioner with ILU subdomain solver
    - 5-level overlap between subdomains
    - Maximum 100 iterations with relative tolerance of 1e-15

    Examples
    --------
    >>> opts = iluasm_options()
    >>> opts = iluasm_options(ksp_max_it=200)
    >>> opts = iluasm_options(pc_asm_overlap=10, ksp_rtol=1e-12)
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
    """
    Set the ISSM MUMPS solver options.

    This function configures MUMPS (MUltifrontal Massively Parallel sparse direct Solver)
    options for use with the ISSM toolkit. It provides default settings optimized for
    ISSM's native solver interface and allows user customization through keyword arguments.

    Parameters
    ----------
    **kwargs : dict
        Additional MUMPS options to override or supplement the defaults.

    Returns
    -------
    collections.OrderedDict
        Dictionary containing ISSM MUMPS solver configuration options with defaults
        updated with any user-provided options.

    Notes
    -----
    This function sets up ISSM's native MUMPS solver interface with:
    - 'mpiparse' matrix type for ISSM's sparse matrix format
    - 'mpi' vector type for parallel vector operations
    - 'mumps' solver type for direct factorization

    Examples
    --------
    >>> opts = issm_mumps_solver()
    >>> opts = issm_mumps_solver(solver_type='mumps')
    """

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
    """
    Set the GSL solver options for ISSM linear solver configuration.

    This function configures GSL (GNU Scientific Library) solver options for use with
    the ISSM toolkit. It provides default settings optimized for ISSM's native GSL
    solver interface and allows user customization through keyword arguments.

    Parameters
    ----------
    **kwargs : dict
        Additional GSL solver options to override or supplement the defaults.

    Returns
    -------
    collections.OrderedDict
        Dictionary containing GSL solver configuration options with defaults
        updated with any user-provided options.

    Notes
    -----
    This function sets up ISSM's native GSL solver interface with:
    - 'dense' matrix type for dense matrix operations
    - 'seq' vector type for sequential vector operations
    - 'gsl' solver type for GSL-based solving

    Examples
    --------
    >>> opts = issm_gsl_solver()
    >>> opts = issm_gsl_solver(solver_type='gsl')
    """

    ## Define defaults
    defaults = {'toolkit': 'issm',
                'mat_type': 'dense',
                'vec_type': 'seq',
                'solver_type': 'gsl'}
    
    ## Update with user options
    gsl = collections.OrderedDict(defaults)
    gsl.update(kwargs)

    return gsl

def is_pc():
    """
    Check if the current operating system is Windows.

    This function determines whether the code is running on a Windows operating system.

    Returns
    -------
    bool
        True if the operating system is Windows, False otherwise.

    Examples
    --------
    >>> is_windows = is_pc()
    >>> print(is_windows)
    True
    """

    return 'Windows' in platform.system()

def get_issm_dir():
    """
    Get the root directory of the ISSM installation.

    This function retrieves the root directory path of the Ice Sheet System Model (ISSM)
    installation by checking the appropriate environment variable based on the operating
    system (ISSM_DIR_WIN for Windows, ISSM_DIR for other systems).

    Returns
    -------
    str
        The root directory path of the ISSM installation.

    Raises
    ------
    EnvironmentError
        If the ISSM_DIR environment variable is not set or is empty.

    Examples
    --------
    >>> issm_dir = get_issm_dir()
    >>> print(issm_dir)
    '/usr/local/issm'
    """

    # Determine ISSM_DIR based on operating system
    if is_pc():
        issm_dir = os.environ['ISSM_DIR_WIN']
        if issm_dir.endswith('/') or issm_dir.endswith('\\'):
            issm_dir = issm_dir[:-1]
    else:
        issm_dir = os.environ['ISSM_DIR']

    # Raise error if ISSM_DIR is not set
    if not issm_dir:
        raise EnvironmentError("ISSM_DIR environment variable is not set. Please ensure ISSM is properly installed and configured.")

    return issm_dir

def get_hostname():
    """
    Get the hostname of the current machine.

    This function retrieves the hostname of the machine on which the code is running.

    Returns
    -------
    str
        The hostname of the current machine.

    Examples
    --------
    >>> hostname = get_hostname()
    """

    return socket.gethostname().lower()

def get_username():
    """
    Get the username of the current user.

    This function retrieves the username of the user currently logged into the system
    by querying the appropriate environment variable based on the operating system.

    Returns
    -------
    str
        The username of the current user. Returns an empty string if the username
        cannot be determined.

    Examples
    --------
    >>> username = get_username()
    >>> print(username)
    'john_doe'
    """

    # Determine username based on operating system
    if 'Windows' in platform.system():
        user_name = os.environ['USERNAME']
    else:
        user_name = os.environ['USER']

    # Return empty string if user_name is None
    if not user_name:
        user_name = ''

    return user_name