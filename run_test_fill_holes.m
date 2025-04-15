
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
        