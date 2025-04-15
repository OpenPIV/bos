import numpy as np

def Parameters(*args):
    """
    Parameters for the BOS (Background-Oriented Schlieren) application
    
    Returns:
    --------
    Mconversion : float
        Conversion factor from pixels to cm
    Const : float
        Constant for the Poisson equation
    Lx : int
        Width of the image in pixels
    Lz : int
        Height of the image in pixels
    val_up : float
        Dirichlet condition at the top
    val_down : float
        Dirichlet condition at the bottom
    nx_pixel : int
        Width of the interrogation window for PIV
    ny_pixel : int
        Height of the interrogation window for PIV
    overlap_x : float
        Overlap ratio in x direction for PIV
    overlap_y : float
        Overlap ratio in y direction for PIV
    """
    # From the calibration body: pixel/cm
    M = 0.0123
    
    # (see the numbers in the BOS sketch)
    B = 0 * (1/M)             # Distance section 5-6  
    W = 20 * (1/M)            # Distance section 3-4
    L = 110 * (1/M)           # Distance section 1-2
    t = 0.5 * (1/M)           # thickness of the glass
    ZD = (L + 2*t + W + B) * (1/M)
    
    # Set the indexes of refraction
    n_air = 1                
    n_water = 1.332
    n_glass = 1.43
    
    # Computing the constant 
    Const_inv = 2 * (((L**2 + B**2)/(n_air)) + (t**2/n_glass) + ((W)/(n_water)))
    Const = (Const_inv)**-1
    
    # Dirichlet's conditions at the top and bottom
    val_up = 1.332    
    val_down = 1.433
    
    # size of the images
    Lx = 1720  # IN PIXEL 
    Lz = 2304
    
    Mconversion = M
    
    # PIV-parameters
    nx_pixel = 64
    ny_pixel = 64
    overlap_x = 0.5
    overlap_y = 0.5
    
    return Mconversion, Const, Lx, Lz, val_up, val_down, nx_pixel, ny_pixel, overlap_x, overlap_y
