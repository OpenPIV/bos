import numpy as np
from openpiv_python import openpiv

def BOS_correlation_OpenPIV(im1, im2, nx_pixel, ny_pixel, overlap_x, overlap_y=None):
    """
    Perform PIV correlation between two images using OpenPIV
    
    Parameters:
    -----------
    im1 : ndarray
        First image
    im2 : ndarray
        Second image
    nx_pixel : int
        Width of the interrogation window
    ny_pixel : int
        Height of the interrogation window
    overlap_x : float
        Overlap ratio in x direction (0-1)
    overlap_y : float, optional
        Overlap ratio in y direction (0-1). If None, uses overlap_x.
        
    Returns:
    --------
    result : dict
        Dictionary containing displacement field data (x, y, u, v)
    """
    # If overlap_y is not provided, use overlap_x
    if overlap_y is None:
        overlap_y = overlap_x
    
    # Calculate overlap in pixels
    ovlapHor = int(nx_pixel * (1 - overlap_x))
    ovlapVer = int(ny_pixel * (1 - overlap_y))
    
    # Call OpenPIV function
    x, y, u, v = openpiv(im1, im2, nx_pixel, ny_pixel, ovlapHor, ovlapVer)
    
    # Create result dictionary
    result = {
        'x': x,
        'y': y,
        'u': u,
        'v': v
    }
    
    return result
