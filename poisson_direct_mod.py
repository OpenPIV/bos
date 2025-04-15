import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve
from scaledata import scaledata

def poisson_direct_mod(Nx, Nz, dxb, dxc, dzb, dzc, rhs):
    """
    Solve the Poisson equation using a direct method with a 5-point stencil
    
    Parameters:
    -----------
    Nx : int
        Number of grid points in x-direction
    Nz : int
        Number of grid points in z-direction
    dxb : ndarray
        Grid spacing in x-direction at boundaries
    dxc : ndarray
        Grid spacing in x-direction at cell centers
    dzb : ndarray
        Grid spacing in z-direction at boundaries
    dzc : ndarray
        Grid spacing in z-direction at cell centers
    rhs : ndarray
        Right-hand side of the Poisson equation
        
    Returns:
    --------
    n2 : ndarray
        Solution of the Poisson equation (refractive index field)
    """
    # Reshape the index array
    idx = np.arange((Nx+2)*(Nz+2)).reshape(Nx+2, Nz+2)
    
    # Initialize arrays for sparse matrix construction
    ii = np.zeros(5*Nx*Nz, dtype=int)
    jj = np.zeros(5*Nx*Nz, dtype=int)
    vv = np.zeros(5*Nx*Nz)
    
    # Compute the terms of the matrix [L] using a stencil of 5 points
    a = 0
    ss = 5
    for k in range(1, Nz+1):
        for i in range(1, Nx+1):
            ii[a:a+ss] = idx[i, k] * np.ones(ss, dtype=int)
            jj[a:a+ss] = np.array([idx[i-1, k], idx[i, k-1], idx[i, k], 
                                   idx[i, k+1], idx[i+1, k]], dtype=int)
            vv[a:a+ss] = np.array([1/(dxc[i]*dxb[i-1]),
                                   1/(dzc[k]*dzb[k-1]),
                                   -1/(dxc[i]*dxb[i-1]) - 1/(dxc[i]*dxb[i]) 
                                   - 1/(dzc[k]*dzb[k-1]) - 1/(dzc[k]*dzb[k]),
                                   1/(dzc[k]*dzb[k]),
                                   1/(dxc[i]*dxb[i])])
            a = a + ss
    
    # Create the sparse matrix L
    L = sparse.csr_matrix((vv, (ii, jj)), shape=((Nx+2)*(Nz+2), (Nx+2)*(Nz+2)))
    
    # Add boundary conditions
    a = 0
    ss = 2
    
    # Bottom boundary
    for i in range(Nx+2):
        ii[a:a+ss] = idx[i, 0] * np.ones(ss, dtype=int)
        jj[a:a+ss] = np.array([idx[i, 0], idx[i, 1]], dtype=int)
        vv[a:a+ss] = np.array([-1, 1])
        a = a + ss
    
    # Top boundary
    for i in range(Nx+2):
        ii[a:a+ss] = idx[i, -1] * np.ones(ss, dtype=int)
        jj[a:a+ss] = np.array([idx[i, -2], idx[i, -1]], dtype=int)
        vv[a:a+ss] = np.array([-1, 1])
        a = a + ss
    
    # Left boundary
    for k in range(1, Nz+1):
        ii[a:a+ss] = idx[0, k] * np.ones(ss, dtype=int)
        jj[a:a+ss] = np.array([idx[0, k], idx[1, k]], dtype=int)
        vv[a:a+ss] = np.array([-1, 1])
        a = a + ss
    
    # Right boundary
    for k in range(1, Nz+1):
        ii[a:a+ss] = idx[-1, k] * np.ones(ss, dtype=int)
        jj[a:a+ss] = np.array([idx[-2, k], idx[-1, k]], dtype=int)
        vv[a:a+ss] = np.array([-1, 1])
        a = a + ss
    
    # Create the boundary condition matrix
    Bp = sparse.csr_matrix((vv[:a], (ii[:a], jj[:a])), shape=((Nx+2)*(Nz+2), (Nx+2)*(Nz+2)))
    
    # Combine the matrices
    L = L + Bp
    
    # Solve the system
    p2 = spsolve(L, rhs.flatten()).reshape(Nx+2, Nz+2)
    
    # Scale the solution
    pg2 = p2 / (np.max(p2) - np.min(p2))
    PG2_gray = pg2 * 255
    n_max = 1.43
    n_min = 1.332
    n2 = scaledata(PG2_gray, n_min, n_max)
    
    return n2
