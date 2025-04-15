import numpy as np
from scipy import signal
from scipy import ndimage
from openpiv_python.cross_correlate_rect import cross_correlate_rect
from openpiv_python.find_displacement_rect import find_displacement_rect
from openpiv_python.sub_pixel_velocity_rect import sub_pixel_velocity_rect
from openpiv_python.fill_holes import fill_holes

def openpiv(a1, b1, ittWidth, ittHeight, ovlapHor, ovlapVer, s2ntype=1, s2nl=1, sclt=1, dt=1, outl=10):
    """
    OpenPIV cross-correlation routine
    
    Parameters:
    -----------
    a1 : ndarray
        First image
    b1 : ndarray
        Second image
    ittWidth : int
        Width of the interrogation window
    ittHeight : int
        Height of the interrogation window
    ovlapHor : int
        Overlap in pixels in horizontal direction
    ovlapVer : int
        Overlap in pixels in vertical direction
    s2ntype : int, optional
        Signal-to-noise ratio type
    s2nl : int, optional
        Signal-to-noise ratio level
    sclt : float, optional
        Scale factor for pixels to physical units
    dt : float, optional
        Time step between images
    outl : float, optional
        Outlier threshold
        
    Returns:
    --------
    x : ndarray
        x-coordinates of the velocity vectors
    y : ndarray
        y-coordinates of the velocity vectors
    u : ndarray
        x-component of the velocity vectors
    v : ndarray
        y-component of the velocity vectors
        
    Notes:
    ------
    How to cite: 
    Taylor, Z.J.; Gurka, R.; Kopp, G.A.; Liberzon, A.; ,
    "Long-Duration Time-Resolved PIV to Study Unsteady Aerodynamics,"
    Instrumentation and Measurement, IEEE Transactions on , vol.59, no.12,
    pp.3262-3269, Dec. 2010 doi: 10.1109/TIM.2010.2047149
    URL: http://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=5464317&isnumber=5609237
    """
    # Ensure images are float32 (equivalent to Matlab's single)
    a1 = a1.astype(np.float32)
    b1 = b1.astype(np.float32)
    
    # Get image dimensions
    verSize, horSize = a1.shape
    
    # Prepare the results storage
    numcols = int(np.floor((horSize - ittWidth) / ovlapHor + 1))
    numrows = int(np.floor((verSize - ittHeight) / ovlapVer + 1))
    res = np.zeros((numcols * numrows, 5))
    resind = 0
    
    # Set origin
    origin = [0, 0]  # origin is bottom left
    
    # FFT size
    NfftWidth = 2 * ittWidth
    NfftHeight = 2 * ittHeight
    
    # Process each interrogation window
    for m in range(0, verSize - ittHeight + 1, ovlapVer):  # vertically
        for k in range(0, horSize - ittWidth + 1, ovlapHor):  # horizontally
            # Extract interrogation windows
            a2 = a1[m:m+ittHeight, k:k+ittWidth]
            b2 = b1[m:m+ittHeight, k:k+ittWidth]
            
            # Cross-correlate
            c = cross_correlate_rect(a2, b2, NfftHeight, NfftWidth)
            
            # Skip if correlation is all zeros
            if not np.any(c):
                u = 0
                v = 0
                y = origin[1] + m + ittHeight/2 - 1
                x = origin[0] + k + ittWidth/2 - 1
                continue
            
            # Find displacement
            peak1, peak2, pixi, pixj = find_displacement_rect(c, s2ntype)
            
            # Sub-pixel interpolation
            peakVer, peakHor, s2n = sub_pixel_velocity_rect(c, pixi, pixj, peak1, peak2, s2nl, ittWidth, ittHeight)
            
            # Scale the pixel displacement to the velocity
            u = (ittWidth - peakHor)
            v = (ittHeight - peakVer)
            y = origin[1] + m + ittHeight/2 - 1
            x = origin[0] + k + ittWidth/2 - 1
            
            # Store results
            resind += 1
            res[resind-1, :] = [x, y, u, v, s2n]
    
    # Reshape U and V matrices in two-dimensional grid
    u = np.reshape(res[:, 2], (numrows, numcols))
    v = np.reshape(res[:, 3], (numrows, numcols))
    vector = u + 1j * v
    
    # Remove outliers - GLOBAL FILTERING
    ind = np.isfinite(np.abs(vector)) & (np.abs(vector) != 0)
    if np.any(ind):
        limit = np.mean(np.abs(vector[ind])) * outl
        outliers = np.abs(vector) > limit
        vector[outliers] = 0
    
    u = np.real(vector)
    v = np.imag(vector)
    
    # Adaptive Local Median filtering
    kernel = np.array([[-1, -1, -1], [-1, 8, -1], [-1, -1, -1]])
    tmpv = np.abs(signal.convolve2d(v, kernel, mode='same'))
    tmpu = np.abs(signal.convolve2d(u, kernel, mode='same'))
    
    # Set limits for outlier detection
    ind = np.isfinite(np.abs(tmpv)) & (np.abs(tmpv) != 0)
    if np.any(ind):
        lmtv = np.mean(tmpv[ind]) + 3 * np.std(tmpv[ind])
    else:
        lmtv = 0
        
    ind = np.isfinite(np.abs(tmpu)) & (np.abs(tmpu) != 0)
    if np.any(ind):
        lmtu = np.mean(tmpu[ind]) + 3 * np.std(tmpu[ind])
    else:
        lmtu = 0
    
    # Find outliers
    u_out = np.where(tmpu > lmtu)
    v_out = np.where(tmpv > lmtv)
    
    # Remove outliers
    u[u_out] = 0
    u[v_out] = 0
    v[v_out] = 0
    v[u_out] = 0
    vector = u + 1j * v
    
    # Update results
    res[:, 2] = np.reshape(np.real(vector), (numrows * numcols))
    res[:, 3] = np.reshape(np.imag(vector), (numrows * numcols))
    
    # Fill holes
    vector = fill_holes(vector)
    res[:, 2] = np.reshape(np.real(vector), (numrows * numcols))
    res[:, 3] = np.reshape(np.imag(vector), (numrows * numcols))
    
    # Scale the pixels and apply the dt
    if sclt != 0:
        res = res * sclt  # pixels to meters
    
    if dt != 0:
        res[:, 2:4] = res[:, 2:4] / dt
    
    # Reshape results
    x = np.reshape(res[:, 0], (numrows, numcols)).T
    y = np.reshape(res[:, 1], (numrows, numcols)).T
    u = np.reshape(res[:, 2], (numrows, numcols)).T
    v = np.reshape(res[:, 3], (numrows, numcols)).T
    
    return x, y, u, v
