import numpy as np

def sub_pixel_velocity_rect(c, pixi, pixj, peak1, peak2, s2nl, ittWidth, ittHeight):
    """
    Calculate Signal-To-Noise Ratio, fit Gaussian bell, find sub-pixel displacement

    Parameters:
    -----------
    c : ndarray
        Cross-correlation matrix
    pixi : int
        Row index of the peak
    pixj : int
        Column index of the peak
    peak1 : float
        Highest peak value
    peak2 : float
        Second highest peak value (or mean value)
    s2nl : float
        Signal-to-noise ratio limit
    ittWidth : int
        Width of the interrogation window
    ittHeight : int
        Height of the interrogation window

    Returns:
    --------
    peakx : float
        Sub-pixel row position of the peak
    peaky : float
        Sub-pixel column position of the peak
    s2n : float
        Signal-to-noise ratio

    Notes:
    ------
    Original MATLAB implementation by Alex Liberzon & Roi Gurka
    Date: Jul-20-99
    """
    # If peak2 equals to zero, it means that nothing was found,
    # and we'll divide by zero
    if peak2 == 0:
        s2n = float('inf')  # Just to protect from zero dividing
    else:
        s2n = peak1 / peak2

    # If Signal-To-Noise ratio is lower than the limit, "mark" it
    if s2n < s2nl:
        peakx = ittHeight
        peaky = ittWidth
    else:  # Otherwise, calculate the velocity
        # Sub-pixel displacement definition by means of Gaussian bell
        if pixi < 2 or pixi > c.shape[0] - 3 or pixj < 2 or pixj > c.shape[1] - 3:
            peakx = ittHeight
            peaky = ittWidth
            return peakx, peaky, s2n

        try:
            # Fit Gaussian in x-direction
            f0 = np.log(max(1e-10, c[pixi, pixj]))
            f1 = np.log(max(1e-10, c[pixi-1, pixj]))
            f2 = np.log(max(1e-10, c[pixi+1, pixj]))
            denom = 2 * f1 - 4 * f0 + 2 * f2
            if abs(denom) > 1e-10:
                peakx = pixi + (f1 - f2) / denom
            else:
                peakx = pixi

            # Fit Gaussian in y-direction
            f0 = np.log(max(1e-10, c[pixi, pixj]))
            f1 = np.log(max(1e-10, c[pixi, pixj-1]))
            f2 = np.log(max(1e-10, c[pixi, pixj+1]))
            denom = 2 * f1 - 4 * f0 + 2 * f2
            if abs(denom) > 1e-10:
                peaky = pixj + (f1 - f2) / denom
            else:
                peaky = pixj
        except:
            peakx = ittHeight
            peaky = ittWidth

        # Check if the result is real (not complex)
        if not np.isreal(peakx) or not np.isreal(peaky):
            peakx = ittHeight
            peaky = ittWidth

    return peakx, peaky, s2n
