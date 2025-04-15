import numpy as np

def scaledata(datain, minval, maxval):
    """
    Scale the values of a matrix from a user specified minimum to a user specified maximum

    Parameters:
    -----------
    datain : ndarray
        Input data to be scaled
    minval : float
        Minimum value for the output data
    maxval : float
        Maximum value for the output data

    Returns:
    --------
    dataout : ndarray
        Scaled data

    Example:
    --------
    >>> a = np.array([1, 2, 3, 4, 5])
    >>> a_out = scaledata(a, 0, 1)
    >>> print(a_out)
    [0.   0.25 0.5  0.75 1.  ]

    Original author: Aniruddha Kembhavi, July 11, 2007
    """
    # Get min and max values
    datamax = np.max(datain)
    datamin = np.min(datain)

    # Scale data exactly as in MATLAB version
    dataout = (datain - datamin) / (datamax - datamin) * (maxval - minval) + minval

    return dataout
