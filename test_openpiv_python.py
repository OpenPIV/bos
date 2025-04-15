import numpy as np
import matplotlib.pyplot as plt
import os
import scipy.io as sio
from openpiv_python.cross_correlate_rect import cross_correlate_rect
from openpiv_python.fill_holes import fill_holes
from openpiv_python.find_displacement_rect import find_displacement_rect
from openpiv_python.inpaint_nans_simple import inpaint_nans
from openpiv_python.sub_pixel_velocity_rect import sub_pixel_velocity_rect

def test_cross_correlate_rect():
    """Test the Python implementation of cross_correlate_rect against the MATLAB implementation"""
    print("Testing cross_correlate_rect function...")

    # Create test data
    size = 32
    a2 = np.random.rand(size, size)
    b2 = np.random.rand(size, size)
    NfftHeight = 64
    NfftWidth = 64

    # Run Python implementation
    result_py = cross_correlate_rect(a2, b2, NfftHeight, NfftWidth)

    # Save test data for MATLAB
    sio.savemat('test_cross_correlate_rect_data.mat', {'a2': a2, 'b2': b2, 'NfftHeight': NfftHeight, 'NfftWidth': NfftWidth})

    # Create MATLAB script to run the test
    with open('run_test_cross_correlate_rect.m', 'w') as f:
        f.write("""
% Add openpiv directory to path
addpath('./openpiv');

% Load test data
load('test_cross_correlate_rect_data.mat');

% Run MATLAB implementation
result_mat = cross_correlate_rect(a2, b2, NfftHeight, NfftWidth);

% Save results
save('test_cross_correlate_rect_results.mat', 'result_mat');

% Exit MATLAB
exit;
        """)

    # Run MATLAB script
    print("Running MATLAB implementation...")
    os.system('matlab -nodisplay -nosplash -nodesktop -r "run(\'run_test_cross_correlate_rect.m\'); exit;"')

    # Load MATLAB results
    try:
        mat_results = sio.loadmat('test_cross_correlate_rect_results.mat')
        result_mat = mat_results['result_mat']

        # Compare results
        diff = np.abs(result_py - result_mat)
        max_diff = np.max(diff)

        print(f"Maximum difference: {max_diff}")

        # Visualize results
        plt.figure(figsize=(12, 4))

        plt.subplot(131)
        plt.imshow(result_py, cmap='viridis')
        plt.title('Python Result')
        plt.colorbar()

        plt.subplot(132)
        plt.imshow(result_mat, cmap='viridis')
        plt.title('MATLAB Result')
        plt.colorbar()

        plt.subplot(133)
        plt.imshow(diff, cmap='hot')
        plt.title('Difference')
        plt.colorbar()

        plt.tight_layout()
        plt.savefig('cross_correlate_rect_comparison.png')

        if max_diff < 1e-10:
            print("Test PASSED: Python and MATLAB implementations produce the same results!")
        else:
            print("Test FAILED: Python and MATLAB implementations produce different results.")

    except Exception as e:
        print(f"Error loading MATLAB results: {e}")
        print("Make sure MATLAB is installed and accessible from the command line.")

def test_inpaint_nans():
    """Test the Python implementation of inpaint_nans against the MATLAB implementation"""
    print("\nTesting inpaint_nans function...")

    # Create test data
    size = 20
    A = np.random.rand(size, size)

    # Add some NaNs
    nan_mask = np.random.rand(size, size) < 0.2
    A[nan_mask] = np.nan

    # Run Python implementation
    result_py = inpaint_nans(A, method=0)

    # Save test data for MATLAB
    sio.savemat('test_inpaint_nans_data.mat', {'A': A})

    # Create MATLAB script to run the test
    with open('run_test_inpaint_nans.m', 'w') as f:
        f.write("""
% Add openpiv directory to path
addpath('./openpiv');

% Load test data
load('test_inpaint_nans_data.mat');

% Run MATLAB implementation
result_mat = inpaint_nans(A, 0);

% Save results
save('test_inpaint_nans_results.mat', 'result_mat');

% Exit MATLAB
exit;
        """)

    # Run MATLAB script
    print("Running MATLAB implementation...")
    os.system('matlab -nodisplay -nosplash -nodesktop -r "run(\'run_test_inpaint_nans.m\'); exit;"')

    # Load MATLAB results
    try:
        mat_results = sio.loadmat('test_inpaint_nans_results.mat')
        result_mat = mat_results['result_mat']

        # Compare results
        diff = np.abs(result_py - result_mat)
        max_diff = np.max(diff)

        print(f"Maximum difference: {max_diff}")

        # Visualize results
        plt.figure(figsize=(15, 5))

        plt.subplot(141)
        plt.imshow(A, cmap='viridis')
        plt.title('Original with NaNs')
        plt.colorbar()

        plt.subplot(142)
        plt.imshow(result_py, cmap='viridis')
        plt.title('Python Result')
        plt.colorbar()

        plt.subplot(143)
        plt.imshow(result_mat, cmap='viridis')
        plt.title('MATLAB Result')
        plt.colorbar()

        plt.subplot(144)
        plt.imshow(diff, cmap='hot')
        plt.title('Difference')
        plt.colorbar()

        plt.tight_layout()
        plt.savefig('inpaint_nans_comparison.png')

        if max_diff < 1e-10:
            print("Test PASSED: Python and MATLAB implementations produce the same results!")
        else:
            print("Test FAILED: Python and MATLAB implementations produce different results.")

    except Exception as e:
        print(f"Error loading MATLAB results: {e}")
        print("Make sure MATLAB is installed and accessible from the command line.")

def test_find_displacement_rect():
    """Test the Python implementation of find_displacement_rect against the MATLAB implementation"""
    print("\nTesting find_displacement_rect function...")

    # Create test data
    size = 64
    c = np.random.rand(size, size)

    # Create a peak
    peak_i, peak_j = 32, 32
    c[peak_i-2:peak_i+3, peak_j-2:peak_j+3] = 0
    c[peak_i, peak_j] = 10

    s2ntype = 1

    # Run Python implementation
    peak1_py, peak2_py, pixi_py, pixj_py = find_displacement_rect(c, s2ntype)

    # Save test data for MATLAB
    sio.savemat('test_find_displacement_rect_data.mat', {'c': c, 's2ntype': s2ntype})

    # Create MATLAB script to run the test
    with open('run_test_find_displacement_rect.m', 'w') as f:
        f.write("""
% Add openpiv directory to path
addpath('./openpiv');

% Load test data
load('test_find_displacement_rect_data.mat');

% Run MATLAB implementation
[peak1_mat, peak2_mat, pixi_mat, pixj_mat] = find_displacement_rect(c, s2ntype);

% Save results
save('test_find_displacement_rect_results.mat', 'peak1_mat', 'peak2_mat', 'pixi_mat', 'pixj_mat');

% Exit MATLAB
exit;
        """)

    # Run MATLAB script
    print("Running MATLAB implementation...")
    os.system('matlab -nodisplay -nosplash -nodesktop -r "run(\'run_test_find_displacement_rect.m\'); exit;"')

    # Load MATLAB results
    try:
        mat_results = sio.loadmat('test_find_displacement_rect_results.mat')
        peak1_mat = mat_results['peak1_mat'][0][0]
        peak2_mat = mat_results['peak2_mat'][0][0]
        pixi_mat = mat_results['pixi_mat'][0][0]
        pixj_mat = mat_results['pixj_mat'][0][0]

        # Compare results
        print(f"Python peak1: {peak1_py}, MATLAB peak1: {peak1_mat}")
        print(f"Python peak2: {peak2_py}, MATLAB peak2: {peak2_mat}")
        print(f"Python pixi: {pixi_py}, MATLAB pixi: {pixi_mat}")
        print(f"Python pixj: {pixj_py}, MATLAB pixj: {pixj_mat}")

        if (abs(peak1_py - peak1_mat) < 1e-10 and
            abs(peak2_py - peak2_mat) < 1e-10 and
            pixi_py == pixi_mat and
            pixj_py == pixj_mat):
            print("Test PASSED: Python and MATLAB implementations produce the same results!")
        else:
            print("Test FAILED: Python and MATLAB implementations produce different results.")

    except Exception as e:
        print(f"Error loading MATLAB results: {e}")
        print("Make sure MATLAB is installed and accessible from the command line.")

def test_sub_pixel_velocity_rect():
    """Test the Python implementation of sub_pixel_velocity_rect against the MATLAB implementation"""
    print("\nTesting sub_pixel_velocity_rect function...")

    # Create test data
    size = 64
    c = np.random.rand(size, size)

    # Create a peak
    peak_i, peak_j = 32, 32
    c[peak_i-2:peak_i+3, peak_j-2:peak_j+3] = 0
    c[peak_i-1, peak_j] = 2
    c[peak_i, peak_j] = 10
    c[peak_i+1, peak_j] = 2
    c[peak_i, peak_j-1] = 2
    c[peak_i, peak_j+1] = 2

    pixi = peak_i
    pixj = peak_j
    peak1 = 10
    peak2 = 2
    s2nl = 1.5
    ittWidth = 32
    ittHeight = 32

    # Run Python implementation
    peakx_py, peaky_py, s2n_py = sub_pixel_velocity_rect(c, pixi, pixj, peak1, peak2, s2nl, ittWidth, ittHeight)

    # Save test data for MATLAB
    sio.savemat('test_sub_pixel_velocity_rect_data.mat',
                {'c': c, 'pixi': pixi, 'pixj': pixj, 'peak1': peak1, 'peak2': peak2,
                 's2nl': s2nl, 'ittWidth': ittWidth, 'ittHeight': ittHeight})

    # Create MATLAB script to run the test
    with open('run_test_sub_pixel_velocity_rect.m', 'w') as f:
        f.write("""
% Add openpiv directory to path
addpath('./openpiv');

% Load test data
load('test_sub_pixel_velocity_rect_data.mat');

% Run MATLAB implementation
[peakx_mat, peaky_mat, s2n_mat] = sub_pixel_velocity_rect(c, pixi, pixj, peak1, peak2, s2nl, ittWidth, ittHeight);

% Save results
save('test_sub_pixel_velocity_rect_results.mat', 'peakx_mat', 'peaky_mat', 's2n_mat');

% Exit MATLAB
exit;
        """)

    # Run MATLAB script
    print("Running MATLAB implementation...")
    os.system('matlab -nodisplay -nosplash -nodesktop -r "run(\'run_test_sub_pixel_velocity_rect.m\'); exit;"')

    # Load MATLAB results
    try:
        mat_results = sio.loadmat('test_sub_pixel_velocity_rect_results.mat')
        peakx_mat = mat_results['peakx_mat'][0][0]
        peaky_mat = mat_results['peaky_mat'][0][0]
        s2n_mat = mat_results['s2n_mat'][0][0]

        # Compare results
        print(f"Python peakx: {peakx_py}, MATLAB peakx: {peakx_mat}")
        print(f"Python peaky: {peaky_py}, MATLAB peaky: {peaky_mat}")
        print(f"Python s2n: {s2n_py}, MATLAB s2n: {s2n_mat}")

        if (abs(peakx_py - peakx_mat) < 1e-10 and
            abs(peaky_py - peaky_mat) < 1e-10 and
            abs(s2n_py - s2n_mat) < 1e-10):
            print("Test PASSED: Python and MATLAB implementations produce the same results!")
        else:
            print("Test FAILED: Python and MATLAB implementations produce different results.")

    except Exception as e:
        print(f"Error loading MATLAB results: {e}")
        print("Make sure MATLAB is installed and accessible from the command line.")

def test_fill_holes():
    """Test the Python implementation of fill_holes against the MATLAB implementation"""
    print("\nTesting fill_holes function...")

    # Create test data
    size = 20
    vector = np.random.rand(size, size) + 1j * np.random.rand(size, size)

    # Add some zeros
    zero_mask = np.random.rand(size, size) < 0.2
    vector[zero_mask] = 0

    # Run Python implementation
    result_py = fill_holes(vector.copy())

    # Save test data for MATLAB
    sio.savemat('test_fill_holes_data.mat', {'vector': vector, 'reslenx': size, 'resleny': size})

    # Create MATLAB script to run the test
    with open('run_test_fill_holes.m', 'w') as f:
        f.write("""
% Add openpiv directory to path
addpath('./openpiv');

% Load test data
load('test_fill_holes_data.mat');

% Run MATLAB implementation
result_mat = fill_holes(vector, reslenx, resleny);

% Save results
save('test_fill_holes_results.mat', 'result_mat');

% Exit MATLAB
exit;
        """)

    # Run MATLAB script
    print("Running MATLAB implementation...")
    os.system('matlab -nodisplay -nosplash -nodesktop -r "run(\'run_test_fill_holes.m\'); exit;"')

    # Load MATLAB results
    try:
        mat_results = sio.loadmat('test_fill_holes_results.mat')
        result_mat = mat_results['result_mat']

        # Compare results
        diff_real = np.abs(np.real(result_py) - np.real(result_mat))
        diff_imag = np.abs(np.imag(result_py) - np.imag(result_mat))
        max_diff_real = np.max(diff_real)
        max_diff_imag = np.max(diff_imag)

        print(f"Maximum difference (real part): {max_diff_real}")
        print(f"Maximum difference (imaginary part): {max_diff_imag}")

        # Visualize results
        plt.figure(figsize=(15, 10))

        plt.subplot(321)
        plt.imshow(np.abs(vector), cmap='viridis')
        plt.title('Original Vector Magnitude')
        plt.colorbar()

        plt.subplot(322)
        plt.imshow(np.angle(vector), cmap='hsv')
        plt.title('Original Vector Phase')
        plt.colorbar()

        plt.subplot(323)
        plt.imshow(np.abs(result_py), cmap='viridis')
        plt.title('Python Result Magnitude')
        plt.colorbar()

        plt.subplot(324)
        plt.imshow(np.angle(result_py), cmap='hsv')
        plt.title('Python Result Phase')
        plt.colorbar()

        plt.subplot(325)
        plt.imshow(np.abs(result_mat), cmap='viridis')
        plt.title('MATLAB Result Magnitude')
        plt.colorbar()

        plt.subplot(326)
        plt.imshow(np.angle(result_mat), cmap='hsv')
        plt.title('MATLAB Result Phase')
        plt.colorbar()

        plt.tight_layout()
        plt.savefig('fill_holes_comparison.png')

        if max_diff_real < 1e-10 and max_diff_imag < 1e-10:
            print("Test PASSED: Python and MATLAB implementations produce the same results!")
        else:
            print("Test FAILED: Python and MATLAB implementations produce different results.")

    except Exception as e:
        print(f"Error loading MATLAB results: {e}")
        print("Make sure MATLAB is installed and accessible from the command line.")

if __name__ == "__main__":
    test_cross_correlate_rect()
    test_inpaint_nans()
    test_find_displacement_rect()
    test_sub_pixel_velocity_rect()
    test_fill_holes()
