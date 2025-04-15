import numpy as np
from PIL import Image
import os
import tkinter as tk
from tkinter import filedialog

def tiffread(filename=None, img_first=None, img_last=None):
    """
    Read TIFF files, including stacks
    
    Parameters:
    -----------
    filename : str, optional
        Path to the TIFF file. If None, a file dialog will open.
    img_first : int, optional
        Index of the first image to read (1-based). If None, all images are read.
    img_last : int, optional
        Index of the last image to read (1-based). If None, all images are read.
        
    Returns:
    --------
    stack : list
        List of dictionaries, each containing image data and metadata
    img_read : int
        Number of images read
        
    Notes:
    ------
    This is a simplified Python version of the tiffread2.m MATLAB function.
    It uses PIL (Python Imaging Library) to read TIFF files.
    
    Original MATLAB implementation by Francois Nedelec, EMBL, Copyright 1999-2007.
    Python translation by Augment Code.
    """
    # If no filename is provided, open a file dialog
    if filename is None:
        root = tk.Tk()
        root.withdraw()
        filename = filedialog.askopenfilename(
            title="Select image file",
            filetypes=[("Image files", "*.tif;*.stk;*.lsm")]
        )
        if not filename:
            return [], 0
    
    # Set default values for img_first and img_last
    if img_first is None:
        img_first = 1
    if img_last is None:
        img_last = 10000
    
    # Open the TIFF file
    try:
        img = Image.open(filename)
    except:
        # Try with .stk extension if .tif fails
        if filename.endswith('.tif'):
            try:
                filename = filename.replace('.tif', '.stk')
                img = Image.open(filename)
            except:
                raise ValueError(f"File {filename} not found or not a valid TIFF file.")
        else:
            raise ValueError(f"File {filename} not found or not a valid TIFF file.")
    
    # Get the number of frames in the TIFF stack
    n_frames = 1
    try:
        while True:
            img.seek(n_frames)
            n_frames += 1
    except EOFError:
        pass
    
    # Adjust img_last if it's greater than the number of frames
    if img_last > n_frames:
        img_last = n_frames
    
    # Read the images
    stack = []
    img_read = 0
    
    for i in range(img_first-1, img_last):
        try:
            img.seek(i)
            
            # Create a dictionary to store image data and metadata
            img_dict = {
                'filename': os.path.abspath(filename),
                'width': img.width,
                'height': img.height,
                'bits': 8 * len(img.getbands()),
                'info': img.info
            }
            
            # Convert image to numpy array
            img_array = np.array(img)
            
            # Handle different image types
            if len(img.getbands()) == 1:
                img_dict['data'] = img_array
            else:
                # For RGB images, split into separate channels
                if len(img.getbands()) >= 3:
                    img_dict['red'] = img_array[:, :, 0]
                    img_dict['green'] = img_array[:, :, 1]
                    img_dict['blue'] = img_array[:, :, 2]
                
                # If there's an alpha channel
                if len(img.getbands()) == 4:
                    img_dict['alpha'] = img_array[:, :, 3]
            
            stack.append(img_dict)
            img_read += 1
            
        except EOFError:
            break
    
    return stack, img_read
