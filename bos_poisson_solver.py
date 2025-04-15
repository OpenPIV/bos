import numpy as np
from crop_field import crop_field
from create_rhs import create_RHS
from create_grid import CreateGrid
from poisson_equation_2d import Poisson_equation_2D

def BOS_PoissonSolver(Displacement_POisson, Const, Lx, Lz):
    """
    Solve the Poisson equation for BOS (Background-Oriented Schlieren) application
    
    Parameters:
    -----------
    Displacement_POisson : dict
        Dictionary containing displacement field data (x, y, u, v)
    Const : float
        Constant for the Poisson equation
    Lx : int
        Width of the image in pixels
    Lz : int
        Height of the image in pixels
        
    Returns:
    --------
    n2 : ndarray
        Solution of the Poisson equation
    xc : ndarray
        x-coordinates of the grid
    zc : ndarray
        z-coordinates of the grid
    """
    # Crop the field because of the remapping algorithm
    Displ = crop_field(Displacement_POisson, Lx, Lz)
    
    # If you decide do not to crop the field, uncomment this line
    # and comment the function crop_field:
    # Displ = Displacement_POisson
    
    # Create the RHS of the POISSON equation by loading the Displ
    Rhs = create_RHS(Displ)
    
    # Rotate the RHS
    rhs = Const * np.fliplr(Rhs)
    Rx, Ry = rhs.shape
    
    # Create the grid
    xc, zc, dxb, dxc, dzc, dzb, Nx, Nz = CreateGrid(Lx, Lz, Rx, Ry)
    
    # Poisson integration
    # n2 = Poisson_directMod(Nx, Nz, dxb, dxc, dzb, dzc, rhs)
    n2 = Poisson_equation_2D(Lx, Lz, Rhs, Const)
    
    return n2, xc, zc
