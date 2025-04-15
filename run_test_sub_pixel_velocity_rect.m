
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
        