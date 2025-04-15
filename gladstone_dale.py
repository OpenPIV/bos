import numpy as np

def Gladstone_Dale(n2, xc, zc):
    """
    Apply Gladstone-Dale relation to convert refractive index to density
    
    Parameters:
    -----------
    n2 : ndarray
        Refractive index field
    xc : ndarray
        x-coordinates of the grid
    zc : ndarray
        z-coordinates of the grid
        
    Returns:
    --------
    Dens : dict
        Dictionary containing density field data (x, z, f)
    Dens_av : ndarray
        Average density profile
        
    Notes:
    ------
    Gladstone-Dale constant for saline-water from:
    https://books.google.co.il/books?id=DJCKI5qQdiAC&pg=PA119&lpg=PA119&dq=gladstone+dale+constant+water
    """
    # Gladstone-Dale constant [g/mL]
    G = 0.335
    
    # Apply Gladstone-Dale relation
    S_out = (n2 - 1) / G
    
    # Create output dictionary
    Dens = {
        'x': xc,
        'z': zc,
        'f': S_out
    }
    
    # Calculate average density profile
    Dens_av = np.mean(S_out, axis=0)
    
    return Dens, Dens_av
