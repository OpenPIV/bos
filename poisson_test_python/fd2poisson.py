import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve
from scipy.linalg import toeplitz

def fd2poisson(Lx, Lz, Nx, Nz, Rhs):
    """
    Numerical approximation to Poisson's equation over the square [0,Lx]x[0,Lz] with
    Dirichlet boundary conditions. Uses a uniform mesh with (Nx+2)x(Nz+2) total
    points (i.e, Nx x Nz interior grid points).

    Parameters:
    -----------
    Lx : float
        Width of the domain
    Lz : float
        Height of the domain
    Nx : int
        Number of interior grid points in x-direction
    Nz : int
        Number of interior grid points in z-direction
    Rhs : ndarray
        Right-hand side of the Poisson equation

    Returns:
    --------
    n2 : ndarray
        Solution of the Poisson equation
    x : ndarray
        x-coordinates of the grid
    y : ndarray
        y-coordinates of the grid
    """
    # Compute the dx and dz
    hx = Lx / (Nx + 1)
    hz = Lz / (Nz + 1)

    # Create mesh, including boundary points
    x, y = np.meshgrid(np.arange(Nx+1), np.arange(Nz+1))

    # Compute u on the boundary from the Dirichlet boundary condition
    ub = np.zeros((Nx, Nz))

    # North and South boundaries
    ub[0, :] = 1.33  # North boundary
    ub[Nz-1, :] = 1.43  # South boundary

    # Convert ub to a vector using column reordering
    ub = (1 / hx**2) * ub.flatten(order='F')

    # Convert f to a vector using column reordering
    f = Rhs.flatten(order='F')

    # Create the D2x and D2y matrices
    z = np.array([-2, 1] + [0] * (Nz - 2))
    D2x = (1 / hx**2) * sparse.kron(sparse.csr_matrix(toeplitz(z, z)), sparse.eye(Nx))
    D2y = (1 / hz**2) * sparse.kron(sparse.eye(Nz), sparse.csr_matrix(toeplitz(z, z)))

    # Solve the system
    u = spsolve(D2x + D2y, f - ub)

    # Convert u from a column vector to a matrix
    u = u.reshape((Nx, Nz), order='F')

    return u, x, y
