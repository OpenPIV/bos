
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
        