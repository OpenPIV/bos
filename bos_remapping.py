import numpy as np
from scipy.interpolate import interp2d
from imwarp import imwarp

def BOS_Remapping(Displ, im1):
    """
    Remap an image according to displacement field
    
    Parameters:
    -----------
    Displ : dict
        Dictionary containing displacement field data (x, y, u, v)
    im1 : ndarray
        Input image to be remapped
        
    Returns:
    --------
    new_image : ndarray
        Remapped image
        
    Notes:
    ------
    Remapping has two steps:
    1. Interpolate the dx, dy fields to each pixel
    2. Warp the image according to the interpolated dx, dy
    """
    # Note that m is vertical, n is horizontal
    m, n = im1.shape
    X, Y = np.meshgrid(np.arange(n), np.arange(m))
    
    # Interpolate displacement fields to each pixel
    v_interp = interp2d(Displ['x'].T, Displ['y'].T, Displ['v'].T, kind='linear')
    u_interp = interp2d(Displ['x'].T, Displ['y'].T, Displ['u'].T, kind='linear')
    
    vi = np.zeros((m, n))
    ui = np.zeros((m, n))
    
    # Apply interpolation to each pixel
    for i in range(m):
        for j in range(n):
            vi[i, j] = v_interp(X[i, j], Y[i, j])
            ui[i, j] = u_interp(X[i, j], Y[i, j])
    
    # Replace NaN values with zeros
    ui[np.isnan(ui)] = 0
    vi[np.isnan(vi)] = 0
    
    # Rotate the image and then apply the remapping
    im1 = np.flipud(im1)
    
    # Apply imwarp and then rotate back the image
    # This ensures the reference is in the same position
    new_image = np.flipud(imwarp(im1, ui, vi, nopad=True))
    
    return new_image
