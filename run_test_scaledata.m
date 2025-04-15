
% Load test data
load('test_scaledata_data.mat');

% Run MATLAB implementation
result_mat = scaledata(datain, minval, maxval);

% Save results
save('test_scaledata_results.mat', 'result_mat');

% Exit MATLAB
exit;
        