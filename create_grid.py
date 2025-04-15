import numpy as np

def CreateGrid(Lx, Lz, Rx, Ry):
    """
    Create a grid for the Poisson solver
    
    Parameters:
    -----------
    Lx : int
        Width of the domain in pixels
    Lz : int
        Height of the domain in pixels
    Rx : int
        Width of the RHS field
    Ry : int
        Height of the RHS field
        
    Returns:
    --------
    xc : ndarray
        x-coordinates of cell centers
    zc : ndarray
        z-coordinates of cell centers
    dxb : ndarray
        Grid spacing in x-direction at boundaries
    dxc : ndarray
        Grid spacing in x-direction at cell centers
    dzc : ndarray
        Grid spacing in z-direction at cell centers
    dzb : ndarray
        Grid spacing in z-direction at boundaries
    Nx : int
        Number of cells in x-direction
    Nz : int
        Number of cells in z-direction
    """
    # Calculate number of cells
    Nx = Rx - 2
    Nz = Ry - 2
    
    # Create grid coordinates
    xi = np.linspace(0, Lx, Nx + 1)
    xb = xi
    dxc = np.ones(len(xi))
    dxb = np.ones(len(xi))
    
    xc = np.linspace(0, Lx, Nx + 2)
    zi = np.linspace(0, Lz, Nz + 1)
    zb = zi
    zc = np.linspace(0, Lz, Nz + 2)
    dzc = np.ones(len(zi))
    dzb = np.ones(len(zi))
    
    return xc, zc, dxb, dxc, dzc, dzb, Nx, Nz
