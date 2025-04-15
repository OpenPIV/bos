
% Add Poisson_test directory to path
addpath('./Poisson_test');

% Run MATLAB implementation
[Mconversion,Const,Lx,Lz,val_up,val_down,nx_pixel,ny_pixel,overlap_x,overlap_y]=Parameters();

% Save results
save('test_parameters_results.mat', 'Mconversion', 'Const', 'Lx', 'Lz', 'val_up', 'val_down', 'nx_pixel', 'ny_pixel', 'overlap_x', 'overlap_y');

% Exit MATLAB
exit;
        