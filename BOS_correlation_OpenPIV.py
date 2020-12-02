from openpiv.pyprocess import extended_search_area_piv as piv
from openpiv.pyprocess import get_coordinates

def BOS_correlation_OpenPIV(im1,im2,nx,ny,overlap_x):
# function [Displ] = BOS_correlation_OpenPIV(im1,im2,nx,ny,overlap_x)
# #  cross-correlation of two images, using OpenPIV (www.openpiv.net)
# #  Inputs:
# #    im1,im2 - images
# #  Outputs:
# #    Displ - displacement field, dx,dy
# returns dictionary of 'x','y','u','v'

    # addpath('./openpiv') 


    #  Overlap in pixels
    overlap_px = int(nx*overlap_x)    #  pix
    overlap_py = int(nx*overlap_x)    #  pix


    #  Compute the cross-correlation using OpenPIV
    #  Note that openpiv also saves .VEC file in the image folder if you need
    #  it later
    #  loadvec([imfile1,'.vec']) 
    
    vel = piv(im1, im2, window_size=nx, overlap=overlap_px)
    x, y = get_coordinates(im1.shape,search_area_size=nx, overlap=overlap_px)
    

    #  not sure it's needed
#     u(isnan(u)) = 0 
#     v(isnan(v)) = 0 
#     u = medfilt2(u, [3 3])   #  size of the window
#     v = medfilt2(v, [3 3]) 

    
#     Displ.x = x 
#     Displ.y = y  
#     Displ.u = u 
#     Displ.v = v  
    
    return {'x':x,'y':y,'u':vel[0],'v':vel[1]}


