import numpy as np
import matplotlib.pyplot as plt

def plotarrow(x, y, u, v, color='b', scale=1):
    """
    Plot an arrow with a head
    
    Parameters:
    -----------
    x : float
        x-coordinate of the arrow start
    y : float
        y-coordinate of the arrow start
    u : float
        x-component of the arrow
    v : float
        y-component of the arrow
    color : str, optional
        Color of the arrow
    scale : float, optional
        Scale factor for the arrow length
        
    Returns:
    --------
    h1 : Line2D
        Line object for the arrow body
    h2 : Line2D
        Line object for the arrow head
    """
    # Parameters for arrow head
    alpha = 0.33  # Size of arrow head relative to the length of the vector
    beta = 0.33   # Width of the base of the arrow head relative to the length
    
    # Scale the arrow
    u = u * scale
    v = v * scale
    
    # Arrow body
    uu = np.array([x, x + u, np.nan])
    vv = np.array([y, y + v, np.nan])
    h1 = plt.plot(uu, vv, color=color)[0]
    
    # Arrow head
    eps = np.finfo(float).eps  # Small value to avoid division by zero
    hu = np.array([x + u - alpha * (u + beta * (v + eps)),
                  x + u,
                  x + u - alpha * (u - beta * (v + eps)),
                  np.nan])
    hv = np.array([y + v - alpha * (v - beta * (u + eps)),
                  y + v,
                  y + v - alpha * (v + beta * (u + eps)),
                  np.nan])
    h2 = plt.plot(hu, hv, color=color)[0]
    
    return h1, h2
