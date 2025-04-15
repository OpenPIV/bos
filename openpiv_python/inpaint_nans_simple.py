import numpy as np
from scipy import interpolate

def inpaint_nans(A, method=0):
    """
    Interpolate NaN values in an array using a simpler approach
    
    Parameters:
    -----------
    A : ndarray
        Array with NaN values to be filled
    method : int, optional
        Not used in this simplified version
        
    Returns:
    --------
    B : ndarray
        Array with NaN values filled in
    """
    # Get array shape
    n, m = A.shape
    
    # Find NaN elements
    nan_mask = np.isnan(A)
    
    # If there are no NaNs, just return the original array
    if not np.any(nan_mask):
        return A.copy()
    
    # Get coordinates of valid values
    x, y = np.meshgrid(np.arange(m), np.arange(n))
    x_valid = x[~nan_mask]
    y_valid = y[~nan_mask]
    values = A[~nan_mask]
    
    # Create interpolator
    if len(values) > 3:  # Need at least 4 points for cubic interpolation
        f = interpolate.griddata(
            (x_valid, y_valid), values, (x, y), 
            method='cubic', fill_value=np.nan
        )
        
        # If there are still NaNs, try linear interpolation
        if np.any(np.isnan(f)):
            f_linear = interpolate.griddata(
                (x_valid, y_valid), values, (x, y),
                method='linear', fill_value=np.nan
            )
            f[np.isnan(f)] = f_linear[np.isnan(f)]
            
            # If there are still NaNs, use nearest neighbor
            if np.any(np.isnan(f)):
                f_nearest = interpolate.griddata(
                    (x_valid, y_valid), values, (x, y),
                    method='nearest', fill_value=np.nan
                )
                f[np.isnan(f)] = f_nearest[np.isnan(f)]
    else:
        # Not enough points for cubic, use linear
        f = interpolate.griddata(
            (x_valid, y_valid), values, (x, y),
            method='linear', fill_value=np.nan
        )
        
        # If there are still NaNs, use nearest neighbor
        if np.any(np.isnan(f)):
            f_nearest = interpolate.griddata(
                (x_valid, y_valid), values, (x, y),
                method='nearest', fill_value=np.nan
            )
            f[np.isnan(f)] = f_nearest[np.isnan(f)]
    
    return f
