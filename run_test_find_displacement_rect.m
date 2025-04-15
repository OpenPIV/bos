
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
        