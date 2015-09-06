function [Displ] = BOS_correlation_OpenPIV(im1,im2,nx,ny,overlap_x)
% cross-correlation of two images, using OpenPIV (www.openpiv.net)
% Inputs:
%   im1,im2 - images
% Outputs:
%   Displ - displacement field, dx,dy

addpath('./openpiv');


% Overlap in pixels
overlap_px = nx*overlap_x;   % pix
overlap_py = nx*overlap_x;   % pix


% Compute the cross-correlation using OpenPIV
% Note that openpiv also saves .VEC file in the image folder if you need
% it later
% loadvec([imfile1,'.vec']);

[x,y,u,v] = openpiv( im1, im2, ...
    nx, ny, ...
    overlap_px, overlap_py);

% not sure it's needed
u(isnan(u)) = 0;
v(isnan(v)) = 0;
u = medfilt2(u, [3 3]);  % size of the window
v = medfilt2(v, [3 3]);


Displ.x = x;
Displ.y = y; 
Displ.u = u;
Displ.v = v; 


