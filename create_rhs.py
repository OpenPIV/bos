import numpy as np

def create_RHS(Displ):
    """
    Create the right-hand side of the Poisson equation
    
    Parameters:
    -----------
    Displ : dict
        Dictionary containing displacement field data (x, y, u, v)
        
    Returns:
    --------
    Rhs : ndarray
        Right-hand side of the Poisson equation
    Nx : int
        Width of the RHS field
    Nz : int
        Height of the RHS field
        
    Notes:
    ------
    This function computes the central difference fields for the displacement
    components and combines them to form the RHS of the Poisson equation.
    """
    # Extract displacement components
    u = Displ['u']
    w = Displ['v']
    x = Displ['x']
    z = Displ['y']
    
    # Get dimensions
    width, height = u.shape
    
    # Preallocate arrays for derivatives
    du = np.zeros((width-2, height-2))
    dw = np.zeros((width-2, height-2))
    
    # Compute central differences
    for k in range(1, width-1):
        for j in range(1, height-1):
            du[k-1, j-1] = (u[k+1, j] - u[k-1, j]) / 2 * abs(x[1, 1] - x[1, 0])
            dw[k-1, j-1] = (w[k, j+1] - w[k, j-1]) / 2 * abs(z[0, 0] - z[0, 1])
    
    # Combine derivatives to form RHS
    Rhs = du + dw
    
    # Replace NaN values with zeros
    Rhs[np.isnan(Rhs)] = 0
    
    # Get dimensions of RHS
    Nz, Nx = Rhs.shape
    
    return Rhs, Nx, Nz
