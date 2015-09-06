function [n2, xc, zc]=BOS_PoissonSolver(Displacement_POisson,Const,Lx,Lz)

% Crop the field because of the remapping algorithm
[Displ]=crop_field(Displacement_POisson,Lx,Lz);


% If you decide do not to crop the field, uncomment this line
% and comment the function crop_field:

%Displ=Displacement_POisson; 

% Create the RHS of the POISSON equation by loading the Displ
[Rhs]=create_RHS(Displ);

% Rotate the RHS
rhs=Const.*fliplr(Rhs);
[Rx,Ry]=size(rhs);

% Create the grid
[xc zc dxb dxc dzc dzb Nx Nz]= CreateGrid(Lx, Lz,Rx,Ry);

% Poisson integration
%[n2]= Poisson_directMod(Nx, Nz, dxb, dxc, dzb, dzc,rhs);
[n2]=Poisson_equation_2D(Lx,Lz,Rhs,Const);