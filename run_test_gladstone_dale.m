
% Load test data
load('test_gladstone_dale_data.mat');

% Run MATLAB implementation
[Dens, Dens_av] = Gladstone_Dale(n2, xc, zc);

% Save results
save('test_gladstone_dale_results.mat', 'Dens', 'Dens_av');

% Exit MATLAB
exit;
        