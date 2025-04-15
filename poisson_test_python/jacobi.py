import numpy as np

def jacobi(Nx, Nz, Rhs):
    """
    Solve Poisson equation using Jacobi iteration method
    
    Parameters:
    -----------
    Nx : int
        Number of grid points in x-direction
    Nz : int
        Number of grid points in z-direction
    Rhs : ndarray
        Right-hand side of the Poisson equation
        
    Returns:
    --------
    n2 : ndarray
        Solution of the Poisson equation
    """
    # Set tolerance and initialize error
    tol = 1e-6  # tolerance
    err = 1.0   # error
    k = 0       # iteration counter
    
    # Initialize the n-value
    nki = Rhs.copy()
    
    # Iterative Jacobi until convergence
    while err > tol:
        k += 1
        
        # Loop through computational nodes (inside the matrix)
        for j in range(1, Nx):
            for i in range(1, Nz):
                nki[j, i] = 0.25 * (Rhs[j+1, i] + Rhs[j, i+1] + Rhs[j-1, i] + Rhs[j, i-1])
        
        # Calculate the error
        err = np.sqrt(np.sum((nki - Rhs)**2))
        
        # Update n
        Rhs = nki.copy()
    
    return nki
