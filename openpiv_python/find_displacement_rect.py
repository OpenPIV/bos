import numpy as np

def find_displacement_rect(c, s2ntype):
    """
    Find the highest peak in cross-correlation matrix and the second peak (or mean value)
    for signal-to-noise ratio calculation.
    
    Parameters:
    -----------
    c : ndarray
        Cross-correlation matrix
    s2ntype : int
        Method (1 or 2) of S2N ratio calculation
        
    Returns:
    --------
    peak1 : float
        Highest peak value
    peak2 : float
        Second highest peak value (or mean value)
    pixi : int
        Row index of the peak1
    pixj : int
        Column index of the peak1
        
    Notes:
    ------
    Original MATLAB implementation by Alex Liberzon & Roi Gurka
    Date: 20-Jul-99
    """
    # Get dimensions of the correlation matrix
    NfftHeight, NfftWidth = c.shape
    
    # Find the maximum peak
    tmp = np.max(c, axis=0)
    pixj = np.argmax(tmp)
    peak1 = tmp[pixj]
    pixi = np.argmax(c[:, pixj])
    
    # Create temporary matrix without the maximum peak
    tmp = c.copy()
    tmp[pixi, pixj] = 0
    
    # If the peak is found on the border, we should not accept it
    if pixi == 0 or pixj == 0 or pixi == NfftHeight-1 or pixj == NfftWidth-1:
        peak2 = peak1  # We'll not accept this peak later, by means of SNR
    else:
        # Look for the Signal-To-Noise ratio by
        # 1. Peak detectability method: First-to-second peak ratio
        # 2. Peak-to-mean ratio - Signal-to-noise estimation
        
        if s2ntype == 1:  # First-to-second peak ratio
            # Remove 3x3 pixels neighborhood around the peak
            i_min = max(0, pixi-1)
            i_max = min(NfftHeight, pixi+2)
            j_min = max(0, pixj-1)
            j_max = min(NfftWidth, pixj+2)
            tmp[i_min:i_max, j_min:j_max] = np.nan
            
            # Look for the second highest peak
            peak2 = np.nanmax(tmp)
            x2, y2 = np.where(tmp == peak2)
            
            # Only if second peak is within the borders
            if len(x2) > 0 and x2[0] > 0 and y2[0] > 0 and x2[0] < NfftHeight-1 and y2[0] < NfftWidth-1:
                # Look for the clear (global) peak, not for a local maximum
                x2, y2 = x2[0], y2[0]
                i_min = max(0, x2-1)
                i_max = min(NfftHeight, x2+2)
                j_min = max(0, y2-1)
                j_max = min(NfftWidth, y2+2)
                
                while peak2 < np.nanmax(c[i_min:i_max, j_min:j_max]):
                    tmp[x2, y2] = np.nan
                    peak2 = np.nanmax(tmp)
                    x2, y2 = np.where(tmp == peak2)
                    if len(x2) == 0:
                        peak2 = peak1  # Will throw this one out later
                        break
                    x2, y2 = x2[0], y2[0]
                    if x2 == 0 or y2 == 0 or x2 == NfftHeight-1 or y2 == NfftWidth-1:
                        peak2 = peak1  # Will throw this one out later
                        break
                    i_min = max(0, x2-1)
                    i_max = min(NfftHeight, x2+2)
                    j_min = max(0, y2-1)
                    j_max = min(NfftWidth, y2+2)
            else:  # Second peak on the border means "second peak doesn't exist"
                peak2 = peak1
        
        elif s2ntype == 2:  # PEAK-TO-MEAN VALUE RATIO
            peak2 = np.nanmean(np.abs(tmp))
    
    return peak1, peak2, pixi, pixj
