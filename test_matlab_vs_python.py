"""
Test script to compare the results of the Matlab and Python implementations of the BOS code.

This script:
1. Runs the Python implementation of the BOS code
2. Runs the Matlab implementation of the BOS code
3. Compares the results of both implementations

Requirements:
- Python 3.6+
- NumPy
- SciPy
- Matplotlib
- Pillow (PIL)
- MATLAB (accessible from the command line)
"""

import numpy as np
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
from PIL import Image
import os
import scipy.io as sio
import time

# Import Python implementations
from parameters import Parameters
from bos_correlation_openpiv import BOS_correlation_OpenPIV
from bos_remapping import BOS_Remapping
from bos_poisson_solver import BOS_PoissonSolver
from gladstone_dale import Gladstone_Dale

def run_python_implementation():
    """Run the Python implementation of the BOS code"""
    print("Running Python implementation...")
    start_time = time.time()
    
    # LOAD THE IMAGES
    im1 = np.array(Image.open('Data/Air_ref.tif')).astype(float) / 255.0
    im2 = np.array(Image.open('Data/Water_ref.tif')).astype(float) / 255.0
    im3 = np.array(Image.open('Data/4layers.tif')).astype(float) / 255.0
    
    # LOAD THE PARAMETERS FILE
    Mconversion, Const, Lx, Lz, val_up, val_down, nx_pixel, ny_pixel, overlap_x, overlap_y = Parameters()
    
    # CREATE THE CALIBRATION FIELD: Correlation air-water
    Calibration = BOS_correlation_OpenPIV(im1, im2, nx_pixel, ny_pixel, overlap_x)
    
    # APPLY THE REMAPPING
    im3_remapped = BOS_Remapping(Calibration, im3)
    Image.fromarray((im3_remapped * 255).astype(np.uint8)).save('Remapped_4layers_python.tif')
    
    # CORRELATION REFERENCE-REMAPPED
    nx_pixel = 32
    ny_pixel = 32
    overlap_x = 0.25
    Displacement_POisson = BOS_correlation_OpenPIV(im1, im3_remapped, nx_pixel, ny_pixel, overlap_x)
    
    # For comparison between corrected and not-corrected case
    Displ_notcorr = BOS_correlation_OpenPIV(im2, im3, nx_pixel, ny_pixel, overlap_x)
    
    # POISSON INTEGRATION
    n2, xc, zc = BOS_PoissonSolver(Displacement_POisson, Const, Lx, Lz)
    n2_nc, x_nc, z_nc = BOS_PoissonSolver(Displ_notcorr, Const, Lx, Lz)
    
    # Gladstone-Dale conversion
    Dens, Dens_av = Gladstone_Dale(n2, xc, zc)
    Dens_2, Dens_av_2 = Gladstone_Dale(n2_nc, x_nc, z_nc)
    
    # Save results for comparison
    np.savez('python_results.npz',
             Calibration=Calibration,
             im3_remapped=im3_remapped,
             Displacement_POisson=Displacement_POisson,
             Displ_notcorr=Displ_notcorr,
             n2=n2, xc=xc, zc=zc,
             n2_nc=n2_nc, x_nc=x_nc, z_nc=z_nc,
             Dens=Dens, Dens_av=Dens_av,
             Dens_2=Dens_2, Dens_av_2=Dens_av_2)
    
    end_time = time.time()
    print(f"Python implementation completed in {end_time - start_time:.2f} seconds")
    
    return {
        'Calibration': Calibration,
        'im3_remapped': im3_remapped,
        'Displacement_POisson': Displacement_POisson,
        'Displ_notcorr': Displ_notcorr,
        'n2': n2, 'xc': xc, 'zc': zc,
        'n2_nc': n2_nc, 'x_nc': x_nc, 'z_nc': z_nc,
        'Dens': Dens, 'Dens_av': Dens_av,
        'Dens_2': Dens_2, 'Dens_av_2': Dens_av_2
    }

def run_matlab_implementation():
    """Run the Matlab implementation of the BOS code"""
    print("Running Matlab implementation...")
    start_time = time.time()
    
    # Create Matlab script to run the BOS code
    with open('run_matlab_bos.m', 'w') as f:
        f.write("""
% Add paths
addpath('./');
addpath('./openpiv');
addpath('./Poisson_test');

% Run BOS_MAIN_ROUTINE
BOS_MAIN_ROUTINE;

% Save results for comparison
save('matlab_results.mat', 'Calibration', 'im3_remapped', 'Displacement_POisson', ...
     'Displ_notcorr', 'n2', 'xc', 'zc', 'n2_nc', 'x_nc', 'z_nc', ...
     'Dens', 'Dens_av', 'Dens_2', 'Dens_av_2');

% Exit MATLAB
exit;
        """)
    
    # Run Matlab script
    os.system('matlab -nodisplay -nosplash -nodesktop -r "run(\'run_matlab_bos.m\'); exit;"')
    
    # Load Matlab results
    try:
        matlab_results = sio.loadmat('matlab_results.mat')
        end_time = time.time()
        print(f"Matlab implementation completed in {end_time - start_time:.2f} seconds")
        return matlab_results
    except Exception as e:
        print(f"Error loading Matlab results: {e}")
        return None

def compare_results(python_results, matlab_results):
    """Compare the results of the Python and Matlab implementations"""
    print("Comparing results...")
    
    if matlab_results is None:
        print("No Matlab results to compare with.")
        return
    
    # Compare calibration field
    u_py = python_results['Calibration']['u']
    v_py = python_results['Calibration']['v']
    u_mat = matlab_results['Calibration']['u'][0][0]
    v_mat = matlab_results['Calibration']['v'][0][0]
    
    # Calculate differences
    u_diff = np.abs(u_py - u_mat)
    v_diff = np.abs(v_py - v_mat)
    
    print(f"Calibration field u: max diff = {np.max(u_diff):.6f}, mean diff = {np.mean(u_diff):.6f}")
    print(f"Calibration field v: max diff = {np.max(v_diff):.6f}, mean diff = {np.mean(v_diff):.6f}")
    
    # Compare remapped image
    im3_remapped_py = python_results['im3_remapped']
    im3_remapped_mat = matlab_results['im3_remapped']
    im3_diff = np.abs(im3_remapped_py - im3_remapped_mat)
    
    print(f"Remapped image: max diff = {np.max(im3_diff):.6f}, mean diff = {np.mean(im3_diff):.6f}")
    
    # Compare displacement field
    u_py = python_results['Displacement_POisson']['u']
    v_py = python_results['Displacement_POisson']['v']
    u_mat = matlab_results['Displacement_POisson']['u'][0][0]
    v_mat = matlab_results['Displacement_POisson']['v'][0][0]
    
    u_diff = np.abs(u_py - u_mat)
    v_diff = np.abs(v_py - v_mat)
    
    print(f"Displacement field u: max diff = {np.max(u_diff):.6f}, mean diff = {np.mean(u_diff):.6f}")
    print(f"Displacement field v: max diff = {np.max(v_diff):.6f}, mean diff = {np.mean(v_diff):.6f}")
    
    # Compare Poisson solution
    n2_py = python_results['n2']
    n2_mat = matlab_results['n2']
    n2_diff = np.abs(n2_py - n2_mat)
    
    print(f"Poisson solution: max diff = {np.max(n2_diff):.6f}, mean diff = {np.mean(n2_diff):.6f}")
    
    # Compare density field
    dens_py = python_results['Dens']['f']
    dens_mat = matlab_results['Dens']['f'][0][0]
    dens_diff = np.abs(dens_py - dens_mat)
    
    print(f"Density field: max diff = {np.max(dens_diff):.6f}, mean diff = {np.mean(dens_diff):.6f}")
    
    # Visualize differences
    plt.figure(figsize=(15, 10))
    
    plt.subplot(231)
    plt.imshow(im3_remapped_py, cmap='gray')
    plt.title('Python: Remapped Image')
    plt.colorbar()
    
    plt.subplot(232)
    plt.imshow(im3_remapped_mat, cmap='gray')
    plt.title('Matlab: Remapped Image')
    plt.colorbar()
    
    plt.subplot(233)
    plt.imshow(im3_diff, cmap='hot')
    plt.title('Difference: Remapped Image')
    plt.colorbar()
    
    plt.subplot(234)
    plt.imshow(n2_py, cmap='viridis')
    plt.title('Python: Poisson Solution')
    plt.colorbar()
    
    plt.subplot(235)
    plt.imshow(n2_mat, cmap='viridis')
    plt.title('Matlab: Poisson Solution')
    plt.colorbar()
    
    plt.subplot(236)
    plt.imshow(n2_diff, cmap='hot')
    plt.title('Difference: Poisson Solution')
    plt.colorbar()
    
    plt.tight_layout()
    plt.savefig('matlab_vs_python_comparison.png')
    print("Comparison results saved to 'matlab_vs_python_comparison.png'")

def main():
    """Main function"""
    # Run Python implementation
    python_results = run_python_implementation()
    
    # Run Matlab implementation
    matlab_results = run_matlab_implementation()
    
    # Compare results
    compare_results(python_results, matlab_results)

if __name__ == "__main__":
    main()
