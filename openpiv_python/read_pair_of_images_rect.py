import numpy as np
from PIL import Image
import cv2

def read_pair_of_images_rect(image1, image2, cropvec, ittWidth, ittHeight, ovlapHor, ovlapVer):
    """
    Read two images and crop them according to 'cropvec'
    
    Parameters:
    -----------
    image1, image2 : str
        Image file paths
    cropvec : list or ndarray
        4-element list of [left, top, right, bottom] - each value is a number of lines
        of interrogation areas (ittWidth pixels) which should be removed before the analysis
    ittWidth : int
        Width of the interrogation window in pixels
    ittHeight : int
        Height of the interrogation window in pixels
    ovlapHor : int
        Horizontal overlap in pixels
    ovlapVer : int
        Vertical overlap in pixels
        
    Returns:
    --------
    A, B : ndarray
        Original images
    A1, B1 : ndarray
        Cropped images
    origin : list
        Origin coordinates [left, top] in pixels
        
    Notes:
    ------
    Original MATLAB implementation by Alex Liberzon & Roi Gurka
    Date: 20-Jul-99
    """
    origin = [0, 0]
    
    try:
        # Read images
        A = np.array(Image.open(image1))
        B = np.array(Image.open(image2))
    except:
        # If standard reading fails, try alternative method
        # Note: tiffread2 is not implemented here, would need a Python equivalent
        raise ValueError("Failed to read images. Alternative method not implemented.")
    
    # Convert to grayscale if color images
    if len(A.shape) == 3:
        A = cv2.cvtColor(A, cv2.COLOR_RGB2GRAY)
        B = cv2.cvtColor(B, cv2.COLOR_RGB2GRAY)
    
    # Normalize to [0, 1]
    A = A.astype(float) / 255.0
    B = B.astype(float) / 255.0
    
    # Get image dimensions
    verSizeA, horSizeA = A.shape
    verSizeB, horSizeB = B.shape
    
    # A & B matrices HAVE to be of the same size, we take smallest
    verSize = min(verSizeA, verSizeB)
    horSize = min(horSizeA, horSizeB)
    
    # Crop images if cropvec is provided
    if np.any(cropvec):
        top = max(0, round(cropvec[1] / ittHeight))  # top side of the image
        left = max(0, round(cropvec[0] / ittWidth))  # left side of the image
        bottom = max(0, np.floor((verSize - cropvec[1] - cropvec[3]) / ittHeight))  # bottom side of the image
        right = max(0, np.floor((horSize - cropvec[0] - cropvec[2]) / ittWidth))  # right of the image
        
        # Calculate crop indices
        top_idx = 1 + top * ittHeight
        bottom_idx = ovlapVer * np.floor(verSize / ovlapVer) - bottom * ittHeight
        left_idx = 1 + left * ittWidth
        right_idx = ovlapHor * np.floor(horSize / ovlapHor) - right * ittWidth
        
        # Convert to integer indices (Python uses 0-based indexing)
        top_idx = int(top_idx) - 1
        bottom_idx = int(bottom_idx)
        left_idx = int(left_idx) - 1
        right_idx = int(right_idx)
        
        # Crop images
        A1 = A[top_idx:bottom_idx, left_idx:right_idx]
        B1 = B[top_idx:bottom_idx, left_idx:right_idx]
    else:
        A1 = A
        B1 = B
        left = 0
        top = 0
    
    # Check if cropped images are valid
    if min(A1.shape) < 1:
        raise ValueError('Zero image or too large interrogation windows')
    
    # Set origin
    origin = [left * ittWidth, top * ittHeight]
    
    return A, B, A1, B1, origin
