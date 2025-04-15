"""
BACKGROUND ORIENTED SCHLIEREN APPLIED TO STRATIFIED LIQUID CASES

An extension to background oriented Schlieren (BOS) is proposed in the
following. The extension enables an accurate reconstruction of the
density field in stratified liquid experiments.
The multi-media imaging through air-glass-water-glass-air leads
to an additional aberration that destroys the reconstruction.
A two-step calibration and image remapping transform are the key
components that correct the images through the stratified media and
provide a non-intrusive full-field density measurements of transparent
liquids.
"""

import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
from parameters import Parameters
from bos_correlation_openpiv import BOS_correlation_OpenPIV
from bos_remapping import BOS_Remapping
from bos_poisson_solver import BOS_PoissonSolver
from gladstone_dale import Gladstone_Dale

def main():
    # LOAD THE IMAGES
    # We capture and load three images of the background pattern, through air
    # (im1), water (im2) and a saline stratified solution (im3)
    # (im stands for image).
    
    im1 = np.array(Image.open('Data/Air_ref.tif')).astype(float) / 255.0
    im2 = np.array(Image.open('Data/Water_ref.tif')).astype(float) / 255.0
    im3 = np.array(Image.open('Data/4layers.tif')).astype(float) / 255.0
    
    # LOAD THE PARAMETERS FILE
    # The parameters file contains: the calibration factor (Mconversion)
    # that has to be measured experimentally since it depends by the camera
    # resolution and distance between camera and backgroud dots; 
    # the images size (Lx,Lz) measured in pixels; the boundary coditions
    # (val_up,val_down) involved in the Poisson's integration; the size in
    # pixels (nx_pixel,ny_pixel) of the interrogation area An, used in the
    # PIV cross-correlation algorithm and the overlap of the shifting window 
    # (overlap_x,overlap_y) in the two directions x and y.
    
    Mconversion, Const, Lx, Lz, val_up, val_down, nx_pixel, ny_pixel, overlap_x, overlap_y = Parameters()
    
    # CREATE THE CALIBRATION FIELD: Correlation air-water
    # The calibration is the displacement field $\Delta x$;$\Delta y$
    # obtained correlating the air and water images (im1,im2).
    
    Calibration = BOS_correlation_OpenPIV(im1, im2, nx_pixel, ny_pixel, overlap_x)
    Magn_cal = np.sqrt(Calibration['u']**2 + Calibration['v']**2)
    
    skip = 2  # Skip vectors
    
    # Plot calibration vector field
    plt.figure(figsize=(12, 8))
    
    plt.subplot(121)
    plt.quiver(Calibration['x'][::skip, ::skip], Calibration['y'][::skip, ::skip],
               Calibration['u'][::skip, ::skip], Calibration['v'][::skip, ::skip],
               scale=1)
    plt.axis('equal')
    plt.xlim([np.min(Calibration['x']), np.max(Calibration['x'])])
    plt.ylim([np.min(Calibration['y']), np.max(Calibration['y'])])
    plt.title('Calibration vector field')
    plt.gca().invert_yaxis()
    plt.xlabel('x [px]')
    plt.ylabel('y [px]')
    
    plt.subplot(122)
    plt.contour(Calibration['x'], Calibration['y'], Magn_cal, 50)
    plt.axis('equal')
    plt.xlim([np.min(Calibration['x']), np.max(Calibration['x'])])
    plt.ylim([np.min(Calibration['y']), np.max(Calibration['y'])])
    plt.gca().invert_yaxis()
    plt.title('Calibration magnitude')
    plt.xlabel('x [px]')
    plt.ylabel('y [px]')
    cbar = plt.colorbar()
    cbar.set_label(r'$\sqrt{\Delta x^{2} +\Delta y^{2}} \hspace{0.25cm}[px]$')
    
    plt.tight_layout()
    
    # THE REMAPPING 
    # Background pattern image obtained through the saline stratified solution
    # is remapped using the displacement field which origins are in the optical
    # system and aberrations due to the multi-media (air-glass-water-glass-air)
    # imaging
    
    im3_remapped = BOS_Remapping(Calibration, im3)
    Image.fromarray((im3_remapped * 255).astype(np.uint8)).save('Remapped_4layers.tif')
    
    # CORRELATION REFERENCE-REMAPPED
    # The corrected image (im3_remapped) is correlated with the original 
    # reference image takein in air (im1) and the result is used to construct
    # the Poission equation and to solve it.
    # We suggest to modify size of the An and overlap for a better
    # cross-correlation result since the displacement im1-im3_remapped is one 
    # order of magnitude smaller than im1-im2 displ.
    
    nx_pixel = 32
    ny_pixel = 32
    overlap_x = 0.25
    skip = 5  # Skip vectors
    
    Displacement_POisson = BOS_correlation_OpenPIV(im1, im3_remapped, nx_pixel, ny_pixel, overlap_x)
    
    # Comparison between corrected and not-corrected case
    Displ_notcorr = BOS_correlation_OpenPIV(im2, im3, nx_pixel, ny_pixel, overlap_x)
    
    plt.figure(figsize=(12, 5))
    
    plt.subplot(121)
    plt.quiver(Displ_notcorr['x'][::skip, ::skip], Displ_notcorr['y'][::skip, ::skip],
               Displ_notcorr['u'][::skip, ::skip], Displ_notcorr['v'][::skip, ::skip], 5)
    plt.axis('equal')
    plt.xlim([np.min(Displ_notcorr['x']), np.max(Displ_notcorr['x'])])
    plt.ylim([np.min(Displ_notcorr['y']), np.max(Displ_notcorr['y'])])
    plt.title('Displacement field without correction')
    plt.gca().invert_yaxis()
    
    plt.subplot(122)
    plt.quiver(Displacement_POisson['x'][::skip, ::skip], Displacement_POisson['y'][::skip, ::skip],
               Displacement_POisson['u'][::skip, ::skip], Displacement_POisson['v'][::skip, ::skip], 5)
    plt.axis('equal')
    plt.xlabel('x [px]')
    plt.ylabel('y [px]')
    plt.xlim([np.min(Displacement_POisson['x']), np.max(Displacement_POisson['x'])])
    plt.ylim([np.min(Displacement_POisson['y']), np.max(Displacement_POisson['y'])])
    plt.gca().invert_yaxis()
    plt.title('Displacement field corrected')
    
    plt.tight_layout()
    
    # POISSON INTEGRATION
    # The result of the correlation (im1,im3_remapped) is than integrate
    # through a Poisson's solutor. Eventually by applying the Gladstone-Dale 
    # conversion we have computed the 2D density field.
    
    n2, xc, zc = BOS_PoissonSolver(Displacement_POisson, Const, Lx, Lz)
    n2_nc, x_nc, z_nc = BOS_PoissonSolver(Displ_notcorr, Const, Lx, Lz)
    
    # Gladstone-Dale conversion
    Dens, Dens_av = Gladstone_Dale(n2, xc, zc)
    Dens_2, Dens_av_2 = Gladstone_Dale(n2_nc, x_nc, z_nc)
    
    # GRAPHICAL OUTPUT 
    
    # Comparison Magnitude
    Magnitudo = np.sqrt(Displacement_POisson['u']**2 + Displacement_POisson['v']**2)
    Magnitudo_nc = np.sqrt(Displ_notcorr['u']**2 + Displ_notcorr['v']**2)
    
    plt.figure(figsize=(12, 5))
    
    plt.subplot(121)
    plt.contourf(Displacement_POisson['x'] * Mconversion,
                Displacement_POisson['y'] * Mconversion, Magnitudo, 50)
    plt.axis('equal')
    cbar1 = plt.colorbar()
    plt.xlabel('x [cm]')
    plt.ylabel('y [cm]')
    plt.xlim([np.min(Displacement_POisson['x'] * Mconversion),
             np.max(Displacement_POisson['x'] * Mconversion)])
    plt.ylim([np.min(Displacement_POisson['y'] * Mconversion),
             np.max(Displacement_POisson['y'] * Mconversion)])
    plt.title('Corrected')
    plt.gca().invert_yaxis()
    cbar1.set_label(r'$\sqrt{\Delta x^{2} +\Delta y^{2}} \hspace{0.25cm}[px]$')
    
    plt.subplot(122)
    plt.contourf(Displ_notcorr['x'] * Mconversion, Displ_notcorr['y'] * Mconversion,
                Magnitudo_nc, 50)
    plt.axis('equal')
    cbar2 = plt.colorbar()
    plt.xlabel('x [cm]')
    plt.ylabel('y [cm]')
    plt.xlim([np.min(Displ_notcorr['x'] * Mconversion),
             np.max(Displ_notcorr['x'] * Mconversion)])
    plt.ylim([np.min(Displ_notcorr['y'] * Mconversion),
             np.max(Displ_notcorr['y'] * Mconversion)])
    plt.title('Not-corrected')
    plt.gca().invert_yaxis()
    cbar2.set_label(r'$\sqrt{\Delta x^{2} +\Delta y^{2}} \hspace{0.25cm}[px]$')
    
    plt.tight_layout()
    
    # Results: Corrected Magnitude, Density field, Density profiles
    fig = plt.figure(figsize=(15, 8))
    
    ax1 = fig.add_axes([0.08, 0.35, 0.3, 0.4])
    cont = ax1.contour(Displacement_POisson['x'] * Mconversion,
                      Displacement_POisson['y'] * Mconversion, Magnitudo, 20)
    plt.colorbar(cont, ax=ax1)
    ax1.set_xlabel('x [cm]')
    ax1.set_ylabel('y [cm]')
    ax1.set_xlim([np.min(Displacement_POisson['x'] * Mconversion),
                 np.max(Displacement_POisson['x'] * Mconversion)])
    ax1.set_ylim([np.min(Displacement_POisson['y'] * Mconversion),
                 np.max(Displacement_POisson['y'] * Mconversion)])
    ax1.set_title('Corrected')
    ax1.invert_yaxis()
    
    ax2 = fig.add_axes([0.42, 0.35, 0.3, 0.4])
    pcm = ax2.pcolormesh(Dens['x'] * Mconversion, Dens['z'] * Mconversion, Dens['f'].T, shading='auto')
    plt.colorbar(pcm, ax=ax2)
    ax2.set_xlabel('x [cm]')
    ax2.set_title('Density Corrected')
    ax2.set_aspect('equal')
    ax2.invert_yaxis()
    
    ax3 = fig.add_axes([0.77, 0.4, 0.20, 0.3])
    ax3.plot(Dens_av, Dens['z'] * Mconversion, 'b', linewidth=1.5, label='Corrected')
    ax3.plot(Dens_av_2, Dens_2['z'] * Mconversion, 'b--', linewidth=1.5, label='Not-Corrected')
    ax3.set_xlabel('ρ [g/mL]')
    ax3.set_xlim([0.99, 1.3])
    ax3.set_ylim([0, 28])
    ax3.legend()
    ax3.invert_yaxis()
    
    plt.tight_layout()
    
    # Comparison Poisson solutions
    plt.figure(figsize=(12, 5))
    
    ax1 = plt.subplot(121, projection='3d')
    ax1.plot_surface(xc, zc, n2.T, cmap='viridis')
    ax1.set_xlabel('x [px]')
    ax1.set_ylabel('y [px]')
    ax1.set_zlabel('n')
    ax1.set_title('Corrected')
    
    ax2 = plt.subplot(122, projection='3d')
    ax2.plot_surface(x_nc, z_nc, n2_nc.T, cmap='viridis')
    ax2.set_xlabel('x [px]')
    ax2.set_ylabel('y [px]')
    ax2.set_zlabel('n')
    ax2.set_title('Not-corrected')
    
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()
