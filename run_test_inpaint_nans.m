
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
        