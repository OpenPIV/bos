
% Load test data
load('test_imwarp_data.mat');

% Run MATLAB implementation
result_mat_nopad = imwarp(I, u, v, true);
result_mat_pad = imwarp(I, u, v);

% Save results
save('test_imwarp_results.mat', 'result_mat_nopad', 'result_mat_pad');

% Exit MATLAB
exit;
        