import numpy as np

def crop_field(Displacement_POisson, Lx, Lz):
    """
    Crop the displacement field to remove boundary artifacts from remapping
    
    Parameters:
    -----------
    Displacement_POisson : dict
        Dictionary containing displacement field data (x, y, u, v)
    Lx : int
        Width of the image in pixels
    Lz : int
        Height of the image in pixels
        
    Returns:
    --------
    Displ : dict
        Cropped displacement field
    """
    # Create a copy of the input dictionary
    Displ = Displacement_POisson.copy()
    
    # Adjust y-coordinates
    Minimum = np.min(Displ['y'])
    Displ['y'] = Displ['y'] - abs(Minimum)
    
    # Calculate magnitude of displacement
    Magnitude = np.sqrt(Displ['u']**2 + Displ['v']**2)
    
    # Get dimensions of the field
    a, b = Displ['x'].shape
    
    # Crop the figure because the remapping algorithm
    # creates the external frames (lack in the data)
    Dxx = Lx / a
    Dyy = Lz / b
    nx_pixels_crop = 700  # 250
    ny_pixels_crop = 350  # 200
    Lx = Lx - nx_pixels_crop
    Lz = Lz - ny_pixels_crop
    
    # Number of pixels to crop on each side
    Dx_pixels = round(nx_pixels_crop / Dxx)
    Dy_pixels = round(ny_pixels_crop / Dyy)
    
    # Crop the displacement field
    Displ['x'] = Displ['x'][Dx_pixels:-Dx_pixels, Dy_pixels:-Dy_pixels]
    Displ['y'] = Displ['y'][Dx_pixels:-Dx_pixels, Dy_pixels:-Dy_pixels]
    Displ['u'] = Displ['u'][Dx_pixels:-Dx_pixels, Dy_pixels:-Dy_pixels]
    Displ['v'] = Displ['v'][Dx_pixels:-Dx_pixels, Dy_pixels:-Dy_pixels]
    
    # Calculate cropped magnitude
    Magnitude_crop = np.sqrt(Displ['u']**2 + Displ['v']**2)
    
    return Displ
