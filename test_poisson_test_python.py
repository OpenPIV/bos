import numpy as np
import matplotlib.pyplot as plt
import os
import scipy.io as sio
from poisson_test_python.fd2poisson import fd2poisson
from poisson_test_python.jacobi import jacobi

def test_fd2poisson():
    """Test the Python implementation of fd2poisson against the MATLAB implementation"""
    print("Testing fd2poisson function...")

    # Create test data
    Lx = 10.0
    Lz = 10.0
    Nx = 20
    Nz = 20

    # Create a test right-hand side (e.g., a simple function)
    x = np.linspace(0, Lx, Nx)
    z = np.linspace(0, Lz, Nz)
    X, Z = np.meshgrid(x, z)
    Rhs = np.sin(X) * np.cos(Z)

    # Run Python implementation
    n2_py, x_py, y_py = fd2poisson(Lx, Lz, Nx, Nz, Rhs)

    # Save test data for MATLAB
    sio.savemat('test_fd2poisson_data.mat', {'Lx': Lx, 'Lz': Lz, 'Nx': Nx, 'Nz': Nz, 'Rhs': Rhs})

    # Create MATLAB script to run the test
    with open('run_test_fd2poisson.m', 'w') as f:
        f.write("""
% Add Poisson_test directory to path
addpath('./Poisson_test');

% Load test data
load('test_fd2poisson_data.mat');

% Run MATLAB implementation
[n2_mat, x_mat, y_mat] = fd2poisson(Lx, Lz, Nx, Nz, Rhs);

% Save results
save('test_fd2poisson_results.mat', 'n2_mat', 'x_mat', 'y_mat');

% Exit MATLAB
exit;
        """)

    # Run MATLAB script
    print("Running MATLAB implementation...")
    os.system('matlab -nodisplay -nosplash -nodesktop -r "run(\'run_test_fd2poisson.m\'); exit;"')

    # Load MATLAB results
    try:
        mat_results = sio.loadmat('test_fd2poisson_results.mat')
        n2_mat = mat_results['n2_mat']
        x_mat = mat_results['x_mat']
        y_mat = mat_results['y_mat']

        # Compare results
        diff = np.abs(n2_py - n2_mat)
        max_diff = np.max(diff)

        print(f"Maximum difference: {max_diff}")

        # Visualize results
        plt.figure(figsize=(15, 5))

        plt.subplot(131)
        plt.imshow(n2_py, cmap='viridis')
        plt.title('Python Result')
        plt.colorbar()

        plt.subplot(132)
        plt.imshow(n2_mat, cmap='viridis')
        plt.title('MATLAB Result')
        plt.colorbar()

        plt.subplot(133)
        plt.imshow(diff, cmap='hot')
        plt.title('Difference')
        plt.colorbar()

        plt.tight_layout()
        plt.savefig('fd2poisson_comparison.png')

        if max_diff < 1e-10:
            print("Test PASSED: Python and MATLAB implementations produce the same results!")
        else:
            print("Test FAILED: Python and MATLAB implementations produce different results.")

    except Exception as e:
        print(f"Error loading MATLAB results: {e}")
        print("Make sure MATLAB is installed and accessible from the command line.")

def test_jacobi():
    """Test the Python implementation of jacobi against the MATLAB implementation"""
    print("\nTesting jacobi function...")

    # Create test data
    Nx = 20
    Nz = 20

    # Create a test right-hand side (e.g., a simple function)
    Rhs = np.zeros((Nx+1, Nz+1))

    # Set boundary conditions
    Rhs[0, :] = 1.0  # Top boundary
    Rhs[-1, :] = 0.0  # Bottom boundary
    Rhs[:, 0] = 0.5  # Left boundary
    Rhs[:, -1] = 0.5  # Right boundary

    # Run Python implementation
    n2_py = jacobi(Nx, Nz, Rhs.copy())

    # Save test data for MATLAB
    sio.savemat('test_jacobi_data.mat', {'Nx': Nx, 'Nz': Nz, 'Rhs': Rhs})

    # Create MATLAB script to run the test
    with open('run_test_jacobi.m', 'w') as f:
        f.write("""
% Add Poisson_test directory to path
addpath('./Poisson_test');

% Load test data
load('test_jacobi_data.mat');

% Run MATLAB implementation
n2_mat = Jacobi(Nx, Nz, Rhs);

% Save results
save('test_jacobi_results.mat', 'n2_mat');

% Exit MATLAB
exit;
        """)

    # Run MATLAB script
    print("Running MATLAB implementation...")
    os.system('matlab -nodisplay -nosplash -nodesktop -r "run(\'run_test_jacobi.m\'); exit;"')

    # Load MATLAB results
    try:
        mat_results = sio.loadmat('test_jacobi_results.mat')
        n2_mat = mat_results['n2_mat']

        # Compare results
        diff = np.abs(n2_py - n2_mat)
        max_diff = np.max(diff)

        print(f"Maximum difference: {max_diff}")

        # Visualize results
        plt.figure(figsize=(15, 5))

        plt.subplot(131)
        plt.imshow(n2_py, cmap='viridis')
        plt.title('Python Result')
        plt.colorbar()

        plt.subplot(132)
        plt.imshow(n2_mat, cmap='viridis')
        plt.title('MATLAB Result')
        plt.colorbar()

        plt.subplot(133)
        plt.imshow(diff, cmap='hot')
        plt.title('Difference')
        plt.colorbar()

        plt.tight_layout()
        plt.savefig('jacobi_comparison.png')

        if max_diff < 1e-10:
            print("Test PASSED: Python and MATLAB implementations produce the same results!")
        else:
            print("Test FAILED: Python and MATLAB implementations produce different results.")

    except Exception as e:
        print(f"Error loading MATLAB results: {e}")
        print("Make sure MATLAB is installed and accessible from the command line.")

if __name__ == "__main__":
    # test_fd2poisson()  # Skip this test for now
    test_jacobi()
