import numpy as np
from openpiv_python.inpaint_nans_simple import inpaint_nans

def fill_holes(vector, reslenx=None, resleny=None):
    """
    Fill holes in a vector field

    Parameters:
    -----------
    vector : ndarray
        Vector field with holes (zeros)
    reslenx : int, optional
        Width of the vector field (not used in this implementation)
    resleny : int, optional
        Height of the vector field (not used in this implementation)

    Returns:
    --------
    vector : ndarray
        Vector field with holes filled

    Notes:
    ------
    Uses inpaint_nans to fill holes in the vector field.
    Holes are identified as elements with zero magnitude.

    Original MATLAB implementation by Alex Liberzon (alex.liberzon@gmail.com)
    Aug 11, 2009
    """
    # Replace zeros with NaNs
    vector[np.abs(vector) == 0] = np.nan

    # Fill holes using inpaint_nans
    real_part = inpaint_nans(np.real(vector))
    imag_part = inpaint_nans(np.imag(vector))

    # Combine real and imaginary parts
    vector = real_part + 1j * imag_part

    return vector
