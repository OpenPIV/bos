import numpy as np
from scipy import signal
from scipy import ndimage

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
    vector = fill_holes(vector, numrows, numcols)
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

def cross_correlate_rect(a, b, NfftHeight, NfftWidth):
    """
    Cross-correlate two rectangular interrogation windows
    
    Parameters:
    -----------
    a : ndarray
        First interrogation window
    b : ndarray
        Second interrogation window
    NfftHeight : int
        Height of the FFT
    NfftWidth : int
        Width of the FFT
        
    Returns:
    --------
    c : ndarray
        Cross-correlation matrix
    """
    # Compute FFT of the interrogation windows
    fa = np.fft.fft2(a, (NfftHeight, NfftWidth))
    fb = np.fft.fft2(b, (NfftHeight, NfftWidth))
    
    # Compute cross-correlation
    c = np.real(np.fft.ifft2(fa * np.conj(fb)))
    
    # Shift the result to have the origin at the center
    c = np.fft.fftshift(c)
    
    return c

def find_displacement_rect(c, s2ntype):
    """
    Find the displacement from the correlation matrix
    
    Parameters:
    -----------
    c : ndarray
        Cross-correlation matrix
    s2ntype : int
        Signal-to-noise ratio type
        
    Returns:
    --------
    peak1 : float
        Value of the highest peak
    peak2 : float
        Value of the second highest peak
    pixi : int
        Row index of the highest peak
    pixj : int
        Column index of the highest peak
    """
    # Find the maximum correlation value and its position
    peak1 = np.max(c)
    ind = np.unravel_index(np.argmax(c), c.shape)
    pixi, pixj = ind
    
    # Create a mask to exclude the region around the first peak
    mask = np.ones_like(c)
    mask_size = 3
    i_min = max(0, pixi - mask_size)
    i_max = min(c.shape[0], pixi + mask_size + 1)
    j_min = max(0, pixj - mask_size)
    j_max = min(c.shape[1], pixj + mask_size + 1)
    mask[i_min:i_max, j_min:j_max] = 0
    
    # Find the second highest peak
    c_masked = c * mask
    peak2 = np.max(c_masked)
    
    return peak1, peak2, pixi, pixj

def sub_pixel_velocity_rect(c, pixi, pixj, peak1, peak2, s2nl, ittWidth, ittHeight):
    """
    Calculate sub-pixel displacement using Gaussian peak fit
    
    Parameters:
    -----------
    c : ndarray
        Cross-correlation matrix
    pixi : int
        Row index of the highest peak
    pixj : int
        Column index of the highest peak
    peak1 : float
        Value of the highest peak
    peak2 : float
        Value of the second highest peak
    s2nl : int
        Signal-to-noise ratio level
    ittWidth : int
        Width of the interrogation window
    ittHeight : int
        Height of the interrogation window
        
    Returns:
    --------
    peakVer : float
        Vertical position of the peak with sub-pixel accuracy
    peakHor : float
        Horizontal position of the peak with sub-pixel accuracy
    s2n : float
        Signal-to-noise ratio
    """
    # Calculate signal-to-noise ratio
    s2n = peak1 / peak2
    
    # Get correlation values around the peak
    if pixi > 0 and pixi < c.shape[0]-1:
        f0 = c[pixi-1, pixj]
        f1 = c[pixi, pixj]
        f2 = c[pixi+1, pixj]
        peakVer = pixi + (np.log(f0) - np.log(f2)) / (2 * np.log(f0) - 4 * np.log(f1) + 2 * np.log(f2))
    else:
        peakVer = pixi
    
    if pixj > 0 and pixj < c.shape[1]-1:
        f0 = c[pixi, pixj-1]
        f1 = c[pixi, pixj]
        f2 = c[pixi, pixj+1]
        peakHor = pixj + (np.log(f0) - np.log(f2)) / (2 * np.log(f0) - 4 * np.log(f1) + 2 * np.log(f2))
    else:
        peakHor = pixj
    
    return peakVer, peakHor, s2n

def fill_holes(vector, numrows, numcols):
    """
    Fill holes in the vector field using interpolation
    
    Parameters:
    -----------
    vector : ndarray
        Complex array containing velocity vectors (u + iv)
    numrows : int
        Number of rows in the vector field
    numcols : int
        Number of columns in the vector field
        
    Returns:
    --------
    vector : ndarray
        Filled vector field
    """
    # Reshape vector to 2D
    vector_2d = vector.copy()
    
    # Find zeros (holes)
    zeros = (vector_2d == 0)
    
    # If there are no holes, return the original vector
    if not np.any(zeros):
        return vector
    
    # Create a mask of valid values
    mask = ~zeros
    
    # Get coordinates of valid values
    y_indices, x_indices = np.where(mask)
    
    # Get values at valid coordinates
    values = vector_2d[mask]
    
    # Create a grid for interpolation
    y_grid, x_grid = np.mgrid[0:numrows, 0:numcols]
    
    # Interpolate using nearest neighbor
    filled_vector = np.zeros_like(vector_2d)
    
    # Interpolate real part
    real_values = np.real(values)
    filled_real = ndimage.griddata(
        (y_indices, x_indices), real_values, 
        (y_grid, x_grid), method='nearest'
    )
    
    # Interpolate imaginary part
    imag_values = np.imag(values)
    filled_imag = ndimage.griddata(
        (y_indices, x_indices), imag_values, 
        (y_grid, x_grid), method='nearest'
    )
    
    # Combine real and imaginary parts
    filled_vector = filled_real + 1j * filled_imag
    
    # Only fill holes, keep original valid values
    vector_2d[zeros] = filled_vector[zeros]
    
    return vector_2d
