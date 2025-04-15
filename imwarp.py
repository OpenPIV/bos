import numpy as np

def imwarp(I, u, v, nopad=False):
    """
    Warp image with flow field
    
    Parameters:
    -----------
    I : ndarray
        Input image
    u : ndarray
        Horizontal displacement field
    v : ndarray
        Vertical displacement field
    nopad : bool, optional
        If True, undefined pixels (source pixel outside the boundary) are given by 
        the nearest boundary pixels, otherwise they are set to NaN.
    
    Returns:
    --------
    O : ndarray
        Warped image
    
    Notes:
    ------
    The flow field is to be given in the coordinate system of O, i.e. the operation 
    warps I toward O.
    
    Original Author: Stefan Roth, Department of Computer Science, TU Darmstadt
    Python translation: Augment Code
    """
    # Image size
    sy, sx = I.shape
    
    # Image size w/ padding
    spx = sx + 2
    spy = sy + 2
    
    if nopad:
        # Warped image coordinates
        X, Y = np.meshgrid(np.arange(1, sx+1), np.arange(1, sy+1))
        XI = (X + u).reshape(1, sx * sy)
        YI = (Y + v).reshape(1, sx * sy)
        
        # Bound coordinates to valid region
        XI = np.maximum(1, np.minimum(sx - 1E-6, XI))
        YI = np.maximum(1, np.minimum(sy - 1E-6, YI))
        
        # Perform linear interpolation (faster than scipy.interpolate.interp2d)
        fXI = np.floor(XI).astype(int)
        cXI = np.ceil(XI).astype(int)
        fYI = np.floor(YI).astype(int)
        cYI = np.ceil(YI).astype(int)
        
        alpha_x = XI - fXI
        alpha_y = YI - fYI
        
        # Convert to 0-based indexing for Python
        fXI = fXI - 1
        cXI = cXI - 1
        fYI = fYI - 1
        cYI = cYI - 1
        
        # Linear interpolation
        O = (1 - alpha_x) * (1 - alpha_y) * I[fYI, fXI] + \
            alpha_x * (1 - alpha_y) * I[fYI, cXI] + \
            (1 - alpha_x) * alpha_y * I[cYI, fXI] + \
            alpha_x * alpha_y * I[cYI, cXI]
    else:
        # Pad image with NaNs
        Z = np.pad(I, 1, mode='constant', constant_values=np.nan)
        
        # Warped image coordinates in padded image
        X, Y = np.meshgrid(np.arange(2, sx+2), np.arange(2, sy+2))
        XI = (X + u).reshape(1, sx * sy)
        YI = (Y + v).reshape(1, sx * sy)
        
        # Bound coordinates to valid region
        XI = np.maximum(1, np.minimum(spx - 1E-6, XI))
        YI = np.maximum(1, np.minimum(spy - 1E-6, YI))
        
        # Perform linear interpolation
        fXI = np.floor(XI).astype(int)
        cXI = np.ceil(XI).astype(int)
        fYI = np.floor(YI).astype(int)
        cYI = np.ceil(YI).astype(int)
        
        alpha_x = XI - fXI
        alpha_y = YI - fYI
        
        # Convert to 0-based indexing for Python
        fXI = fXI - 1
        cXI = cXI - 1
        fYI = fYI - 1
        cYI = cYI - 1
        
        # Linear interpolation
        O = (1 - alpha_x) * (1 - alpha_y) * Z[fYI, fXI] + \
            alpha_x * (1 - alpha_y) * Z[fYI, cXI] + \
            (1 - alpha_x) * alpha_y * Z[cYI, fXI] + \
            alpha_x * alpha_y * Z[cYI, cXI]
    
    return O.reshape(sy, sx)
