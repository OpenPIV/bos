import numpy as np

def cross_correlate_rect(a2, b2, NfftHeight, NfftWidth):
    """
    Cross-correlate two rectangular interrogation windows
    
    Parameters:
    -----------
    a2 : ndarray
        First interrogation window
    b2 : ndarray
        Second interrogation window
    NfftHeight : int
        Height of the FFT
    NfftWidth : int
        Width of the FFT
        
    Returns:
    --------
    c : ndarray
        Cross-correlation matrix
        
    Notes:
    ------
    A sort of normalized cross correlation of two images using 
    convolution of conjugates and inverse FFT.
    Only the real and non-negative part is preserved.
    Average intensity of each image is subtracted (a sort of highpass).
    """
    # Subtract mean from each window
    a2 = a2 - np.mean(a2)
    b2 = b2 - np.mean(b2)
    
    # Flip the second window
    b2 = np.flipud(np.fliplr(b2))
    
    # Compute FFT of the windows
    ffta = np.fft.fft2(a2, (NfftHeight, NfftWidth))
    fftb = np.fft.fft2(b2, (NfftHeight, NfftWidth))
    
    # Compute cross-correlation
    c = np.real(np.fft.ifft2(ffta * fftb))
    
    # Set negative values to zero
    c[c < 0] = 0
    
    return c
