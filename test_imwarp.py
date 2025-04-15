import numpy as np
import matplotlib.pyplot as plt
from imwarp import imwarp
import os
import scipy.io as sio

def test_imwarp():
    """
    Test the Python implementation of imwarp against the MATLAB implementation.
    
    This test:
    1. Creates a test image and displacement fields
    2. Saves them as .mat files for MATLAB
    3. Runs the Python implementation
    4. Calls a MATLAB script to run the MATLAB implementation and save results
    5. Loads the MATLAB results and compares with Python results
    """
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
        plt.show()
        
        if max_diff_nopad < 1e-10 and max_diff_pad < 1e-10:
            print("Test PASSED: Python and MATLAB implementations produce the same results!")
        else:
            print("Test FAILED: Python and MATLAB implementations produce different results.")
            
    except Exception as e:
        print(f"Error loading MATLAB results: {e}")
        print("Make sure MATLAB is installed and accessible from the command line.")

if __name__ == "__main__":
    test_imwarp()
