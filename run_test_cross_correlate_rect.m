
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
        