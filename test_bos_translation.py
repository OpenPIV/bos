import numpy as np
import matplotlib.pyplot as plt
from PIL import Image
import os
import scipy.io as sio
from imwarp import imwarp
from parameters import Parameters
from bos_remapping import BOS_Remapping
from gladstone_dale import Gladstone_Dale
from crop_field import crop_field
from create_rhs import create_RHS
from create_grid import CreateGrid
from poisson_equation_2d import Poisson_equation_2D
from scaledata import scaledata

def test_imwarp():
    """Test the Python implementation of imwarp against the MATLAB implementation"""
    print("Testing imwarp function...")

    # Create test data
    size = 50
    x = np.linspace(-5, 5, size)
    y = np.linspace(-5, 5, size)
    X, Y = np.meshgrid(x, y)

    # Create a test image (a Gaussian)
    I = np.exp(-(X**2 + Y**2) / 2)

    # Create displacement fields (a simple rotation-like field)
    u = 0.1 * Y
    v = -0.1 * X

    # Run Python implementation
    result_py_nopad = imwarp(I, u, v, nopad=True)
    result_py_pad = imwarp(I, u, v, nopad=False)

    # Save test data for MATLAB
    sio.savemat('test_imwarp_data.mat', {'I': I, 'u': u, 'v': v})

    # Create MATLAB script to run the test
    with open('run_test_imwarp.m', 'w') as f:
        f.write("""
% Load test data
load('test_imwarp_data.mat');

% Run MATLAB implementation
result_mat_nopad = imwarp(I, u, v, true);
result_mat_pad = imwarp(I, u, v);

% Save results
save('test_imwarp_results.mat', 'result_mat_nopad', 'result_mat_pad');

% Exit MATLAB
exit;
        """)

    # Run MATLAB script
    print("Running MATLAB implementation...")
    os.system('matlab -nodisplay -nosplash -nodesktop -r "run(\'run_test_imwarp.m\'); exit;"')

    # Load MATLAB results
    try:
        mat_results = sio.loadmat('test_imwarp_results.mat')
        result_mat_nopad = mat_results['result_mat_nopad']
        result_mat_pad = mat_results['result_mat_pad']

        # Compare results
        diff_nopad = np.abs(result_py_nopad - result_mat_nopad)
        diff_pad = np.abs(result_py_pad - result_mat_pad)

        max_diff_nopad = np.nanmax(diff_nopad)
        max_diff_pad = np.nanmax(diff_pad)

        print(f"Maximum difference (nopad=True): {max_diff_nopad}")
        print(f"Maximum difference (nopad=False): {max_diff_pad}")

        # Visualize results
        plt.figure(figsize=(15, 10))

        plt.subplot(2, 3, 1)
        plt.imshow(I)
        plt.title('Original Image')

        plt.subplot(2, 3, 2)
        plt.imshow(result_py_nopad)
        plt.title('Python Result (nopad=True)')

        plt.subplot(2, 3, 3)
        plt.imshow(result_mat_nopad)
        plt.title('MATLAB Result (nopad=True)')

        plt.subplot(2, 3, 4)
        plt.quiver(X[::3, ::3], Y[::3, ::3], u[::3, ::3], v[::3, ::3])
        plt.title('Displacement Field')

        plt.subplot(2, 3, 5)
        plt.imshow(result_py_pad)
        plt.title('Python Result (nopad=False)')

        plt.subplot(2, 3, 6)
        plt.imshow(result_mat_pad)
        plt.title('MATLAB Result (nopad=False)')

        plt.tight_layout()
        plt.savefig('imwarp_comparison.png')

        if max_diff_nopad < 1e-10 and max_diff_pad < 1e-10:
            print("Test PASSED: Python and MATLAB implementations produce the same results!")
        else:
            print("Test FAILED: Python and MATLAB implementations produce different results.")

    except Exception as e:
        print(f"Error loading MATLAB results: {e}")
        print("Make sure MATLAB is installed and accessible from the command line.")

def test_parameters():
    """Test the Python implementation of Parameters against the MATLAB implementation"""
    print("\nTesting Parameters function...")

    # Run Python implementation
    Mconversion, Const, Lx, Lz, val_up, val_down, nx_pixel, ny_pixel, overlap_x, overlap_y = Parameters()

    # Create MATLAB script to run the test
    with open('run_test_parameters.m', 'w') as f:
        f.write("""
% Add Poisson_test directory to path
addpath('./Poisson_test');

% Run MATLAB implementation
[Mconversion,Const,Lx,Lz,val_up,val_down,nx_pixel,ny_pixel,overlap_x,overlap_y]=Parameters();

% Save results
save('test_parameters_results.mat', 'Mconversion', 'Const', 'Lx', 'Lz', 'val_up', 'val_down', 'nx_pixel', 'ny_pixel', 'overlap_x', 'overlap_y');

% Exit MATLAB
exit;
        """)

    # Run MATLAB script
    print("Running MATLAB implementation...")
    os.system('matlab -nodisplay -nosplash -nodesktop -r "run(\'run_test_parameters.m\'); exit;"')

    # Load MATLAB results
    try:
        mat_results = sio.loadmat('test_parameters_results.mat')

        # Compare results
        print(f"Python Mconversion: {Mconversion}, MATLAB Mconversion: {mat_results['Mconversion'][0][0]}")
        print(f"Python Const: {Const}, MATLAB Const: {mat_results['Const'][0][0]}")
        print(f"Python Lx: {Lx}, MATLAB Lx: {mat_results['Lx'][0][0]}")
        print(f"Python Lz: {Lz}, MATLAB Lz: {mat_results['Lz'][0][0]}")

        if (abs(Mconversion - mat_results['Mconversion'][0][0]) < 1e-10 and
            abs(Const - mat_results['Const'][0][0]) < 1e-10 and
            Lx == mat_results['Lx'][0][0] and
            Lz == mat_results['Lz'][0][0]):
            print("Test PASSED: Python and MATLAB implementations produce the same results!")
        else:
            print("Test FAILED: Python and MATLAB implementations produce different results.")

    except Exception as e:
        print(f"Error loading MATLAB results: {e}")
        print("Make sure MATLAB is installed and accessible from the command line.")

def test_scaledata():
    """Test the Python implementation of scaledata against the MATLAB implementation"""
    print("\nTesting scaledata function...")

    # Create test data
    datain = np.array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
    minval = 0
    maxval = 1

    # Run Python implementation
    result_py = scaledata(datain, minval, maxval)

    # Save test data for MATLAB
    sio.savemat('test_scaledata_data.mat', {'datain': datain, 'minval': minval, 'maxval': maxval})

    # Create MATLAB script to run the test
    with open('run_test_scaledata.m', 'w') as f:
        f.write("""
% Load test data
load('test_scaledata_data.mat');

% Run MATLAB implementation
result_mat = scaledata(datain, minval, maxval);

% Save results
save('test_scaledata_results.mat', 'result_mat');

% Exit MATLAB
exit;
        """)

    # Run MATLAB script
    print("Running MATLAB implementation...")
    os.system('matlab -nodisplay -nosplash -nodesktop -r "run(\'run_test_scaledata.m\'); exit;"')

    # Load MATLAB results
    try:
        mat_results = sio.loadmat('test_scaledata_results.mat')
        result_mat = mat_results['result_mat']

        # Compare results
        diff = np.abs(result_py - result_mat)
        max_diff = np.max(diff)

        print(f"Python result:\n{result_py}")
        print(f"MATLAB result:\n{result_mat}")
        print(f"Maximum difference: {max_diff}")

        if max_diff < 1e-10:
            print("Test PASSED: Python and MATLAB implementations produce the same results!")
        else:
            print("Test FAILED: Python and MATLAB implementations produce different results.")

    except Exception as e:
        print(f"Error loading MATLAB results: {e}")
        print("Make sure MATLAB is installed and accessible from the command line.")

def test_gladstone_dale():
    """Test the Python implementation of Gladstone_Dale against the MATLAB implementation"""
    print("\nTesting Gladstone_Dale function...")

    # Create test data
    n2 = np.random.rand(10, 10) * 0.1 + 1.3
    xc = np.linspace(0, 1, 10)
    zc = np.linspace(0, 1, 10)

    # Run Python implementation
    Dens_py, Dens_av_py = Gladstone_Dale(n2, xc, zc)

    # Save test data for MATLAB
    sio.savemat('test_gladstone_dale_data.mat', {'n2': n2, 'xc': xc, 'zc': zc})

    # Create MATLAB script to run the test
    with open('run_test_gladstone_dale.m', 'w') as f:
        f.write("""
% Load test data
load('test_gladstone_dale_data.mat');

% Run MATLAB implementation
[Dens, Dens_av] = Gladstone_Dale(n2, xc, zc);

% Save results
save('test_gladstone_dale_results.mat', 'Dens', 'Dens_av');

% Exit MATLAB
exit;
        """)

    # Run MATLAB script
    print("Running MATLAB implementation...")
    os.system('matlab -nodisplay -nosplash -nodesktop -r "run(\'run_test_gladstone_dale.m\'); exit;"')

    # Load MATLAB results
    try:
        mat_results = sio.loadmat('test_gladstone_dale_results.mat')
        Dens_mat = mat_results['Dens']
        Dens_av_mat = mat_results['Dens_av'][0]

        # Compare results
        diff_f = np.abs(Dens_py['f'] - Dens_mat['f'][0][0])
        diff_av = np.abs(Dens_av_py - Dens_av_mat)

        max_diff_f = np.max(diff_f)
        max_diff_av = np.max(diff_av)

        print(f"Maximum difference in density field: {max_diff_f}")
        print(f"Maximum difference in average density: {max_diff_av}")

        if max_diff_f < 1e-10 and max_diff_av < 1e-10:
            print("Test PASSED: Python and MATLAB implementations produce the same results!")
        else:
            print("Test FAILED: Python and MATLAB implementations produce different results.")

    except Exception as e:
        print(f"Error loading MATLAB results: {e}")
        print("Make sure MATLAB is installed and accessible from the command line.")

if __name__ == "__main__":
    test_imwarp()
    test_parameters()
    test_scaledata()
    test_gladstone_dale()
