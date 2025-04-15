import numpy as np
from scaledata import scaledata

def Poisson_equation_2D(Lx, Lz, Rhs, Const):
    """
    Solve the 2D Poisson equation using the Finite Difference Method
    
    Parameters:
    -----------
    Lx : int
        Width of the domain in pixels
    Lz : int
        Height of the domain in pixels
    Rhs : ndarray
        Right-hand side of the Poisson equation
    Const : float
        Constant for the Poisson equation
        
    Returns:
    --------
    n2 : ndarray
        Solution of the Poisson equation (refractive index field)
        
    Notes:
    ------
    Numerical scheme used is a second order central difference in space
    (5-point difference)
    """
    # Get dimensions of RHS
    Nx, Nz = Rhs.shape
    
    # Specify parameters
    dx = Lx / (Nx - 1)  # Width of space step (x)
    dy = Lz / (Nz - 1)  # Width of space step (y)
    x = np.linspace(0, Lx, Nx)  # Range of x and grid points
    y = np.linspace(0, Lz, Nz)  # Range of y and grid points
    
    # Preallocate arrays
    b = np.zeros((Nx, Nz))
    pn = np.zeros((Nx, Nz))
    
    # Initial conditions
    p = np.zeros((Nx, Nz))
    
    # Prepare RHS
    Rhs = Const * np.fliplr(Rhs)
    b = Rhs
    
    # Define interior points
    i = np.arange(1, Nx-1)
    j = np.arange(1, Nz-1)
    
    # Poisson equation solution (iterative method)
    tol = 1e-4  # Set tolerance
    maxerr = float('inf')  # Initial error
    iter_count = 0
    pn = p.copy()
    
    # Iterative solution
    while maxerr > tol:
        iter_count += 1
        
        # Explicit iterative scheme with central difference in space (5-point difference)
        for ii in i:
            for jj in j:
                p[ii, jj] = ((dy**2 * (pn[ii+1, jj] + pn[ii-1, jj])) + 
                             (dx**2 * (pn[ii, jj+1] + pn[ii, jj-1])) - 
                             (b[ii, jj] * dx**2 * dy * 2)) / (2 * (dx**2 + dy**2))
        
        # Boundary conditions
        # Neumann's conditions: dp/dx|end = dp/dx|end-1
        p[0, :] = p[1, :]
        p[-1, :] = p[-2, :]
        
        # Neumann's conditions
        p[:, 0] = p[:, 1]
        p[:, -1] = p[:, -2]
        
        # Calculate error
        with np.errstate(divide='ignore', invalid='ignore'):
            err_matrix = np.abs((p - pn) / p)
            err_matrix[~np.isfinite(err_matrix)] = 0
            maxerr = np.max(err_matrix)
        
        # Update previous solution
        pn = p.copy()
    
    # Scale the solution to physical values
    PG2_gray = p * 255
    n_max = 1.43
    n_min = 1.332
    n2 = scaledata(PG2_gray, n_min, n_max)
    
    return n2
