import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
import os
from parameters import Parameters
from bos_correlation_openpiv import BOS_correlation_OpenPIV
from bos_remapping import BOS_Remapping
from bos_poisson_solver import BOS_PoissonSolver
from gladstone_dale import Gladstone_Dale

def main():
    """
    Main script for the BOS (Background-Oriented Schlieren) application
    
    This script:
    1. Loads reference and test images
    2. Loads parameters
    3. Creates calibration field by correlating air-water reference images
    4. Applies remapping to test image
    5. Correlates reference with remapped image
    6. Compares with non-corrected case
    7. Solves Poisson equation for both cases
    8. Applies Gladstone-Dale conversion
    9. Generates graphical outputs
    """
    # LOAD THE IMAGES
    im1 = np.array(Image.open('Data/Air_ref.tif')).astype(float) / 255.0
    im2 = np.array(Image.open('Data/Water_ref.tif')).astype(float) / 255.0
    im3 = np.array(Image.open('Data/4layers.tif')).astype(float) / 255.0
    
    # LOAD THE PARAMETERS FILE
    Mconversion, Const, Lx, Lz, val_up, val_down, nx_pixel, ny_pixel, overlap_x, overlap_y = Parameters()
    
    # CREATE THE CALIBRATION FIELD: Correlation air-water
    Calibration = BOS_correlation_OpenPIV(im1, im2, nx_pixel, ny_pixel, overlap_x)
    
    # Plot calibration field
    plt.figure()
    plt.quiver(Calibration['x'], Calibration['y'], Calibration['u'], Calibration['v'], scale=1)
    plt.axis('equal')
    
    Magn_cal = np.sqrt(Calibration['u']**2 + Calibration['v']**2)
    plt.figure()
    plt.contour(Calibration['x'], Calibration['y'], Magn_cal, 50)
    plt.axis('equal')
    plt.colorbar()
    
    # APPLY THE REMAPPING
    im3_remapped = BOS_Remapping(Calibration, im3)
    
    # Save remapped image
    Image.fromarray((im3_remapped * 255).astype(np.uint8)).save('Remapped_4layers.tif')
    
    # CORRELATION REFERENCE-REMAPPED
    nx_pixel = 32
    ny_pixel = 32
    overlap_x = 0.25
    Displacement_POisson = BOS_correlation_OpenPIV(im1, im3_remapped, nx_pixel, ny_pixel, overlap_x)
    
    # Check the displacement corrected
    plt.figure()
    plt.quiver(Displacement_POisson['x'], Displacement_POisson['y'], 
               Displacement_POisson['u'], Displacement_POisson['v'], scale=0.2)
    plt.axis('equal')
    
    # For comparison between corrected and not-corrected case
    # Small displacement im2-im3. The interrogation area needs to be reduced
    Displ_notcorr = BOS_correlation_OpenPIV(im2, im3, nx_pixel, ny_pixel, overlap_x)
    
    # Check the displacement not-corrected
    plt.figure()
    plt.quiver(Displ_notcorr['x'], Displ_notcorr['y'], 
               Displ_notcorr['u'], Displ_notcorr['v'], scale=0.2)
    plt.axis('equal')
    
    # POISSON INTEGRATION
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
    plt.contour(Displacement_POisson['x'] * Mconversion, 
                Displacement_POisson['y'] * Mconversion, 
                Magnitudo, 20)
    plt.axis('equal')
    plt.colorbar()
    plt.xlabel('x [cm]')
    plt.ylabel('y [cm]')
    plt.xlim([np.min(Displacement_POisson['x'] * Mconversion), 
              np.max(Displacement_POisson['x'] * Mconversion)])
    plt.ylim([np.min(Displacement_POisson['y'] * Mconversion), 
              np.max(Displacement_POisson['y'] * Mconversion)])
    plt.clim(0, 5)
    plt.title('Corrected')
    plt.gca().invert_yaxis()
    
    plt.subplot(122)
    plt.contour(Displ_notcorr['x'] * Mconversion, 
                Displ_notcorr['y'] * Mconversion, 
                Magnitudo_nc, 20)
    plt.axis('equal')
    plt.colorbar()
    plt.xlabel('x [cm]')
    plt.ylabel('y [cm]')
    plt.xlim([np.min(Displ_notcorr['x'] * Mconversion), 
              np.max(Displ_notcorr['x'] * Mconversion)])
    plt.ylim([np.min(Displ_notcorr['y'] * Mconversion), 
              np.max(Displ_notcorr['y'] * Mconversion)])
    plt.clim(0, 5)
    plt.title('Not-corrected')
    plt.gca().invert_yaxis()
    
    # Results: Corrected Magnitude, Density field, Density profiles
    fig = plt.figure(figsize=(15, 8))
    
    ax1 = fig.add_axes([0.08, 0.35, 0.3, 0.4])
    cont = ax1.contour(Displacement_POisson['x'] * Mconversion, 
                       Displacement_POisson['y'] * Mconversion, 
                       Magnitudo, 20)
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
