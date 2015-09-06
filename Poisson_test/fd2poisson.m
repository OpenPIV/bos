% Numerical approximation to Poisson’s equation over the square [a,b]x[a,b] with
% Dirichlet boundary conditions. Uses a uniform mesh with (n+2)x(n+2) total
% points (i.e, n x n interior grid points).
% Input:
% pfunc : the RHS of poisson equation (i.e. the Laplacian of u).
% bfunc : the boundary function representing the Dirichlet B.C.
% a,b : the interval defining the square
% n : n+2 is the number of points in either direction of the mesh.
% Ouput:
% u : the numerical solution of Poisson equation at the mesh points.
% x,y : the uniform mesh.


function [n2,x,y] = fd2poisson(Lx,Lz,Nx,Nz,Rhs)

% compute the dx and dz
hx = (Lx)/(Nx+1); 
hz = (Lz)/(Nz+1);

[x,y] = meshgrid(0:Nx,0:Nz); % mesh, including boundary points. 

% Compute u on the boundary from the Dirichlet boundary condition
ub = zeros(Nx,Nz);
idx = 2:Nx+1;
idy = 2:Nz+1;

% West and East boundaries need special attention
% ub(:,1) = ; % West Boundary
% ub(:,n) = ; % East Boundary

% Now the North and South boundaries
ub(1,1:Nx) = 1.33;
ub(Nz,1:Nx) = 1.43;
% Convert ub to a vector using column reordering
ub = (1/hx^2)*reshape(ub,Nx*Nz,1);

% Call the RHS of Poisson’s equation at the interior points.
% Convert f to a vector using column reordering
f = reshape(Rhs,Nx*Nz,1);

% Create the D2x and D2y matrices
z = [-2;1;zeros(Nz-2,1)];
D2x = 1/h^2*kron(toeplitz(z,z),eye(Nx));
D2y = 1/h^2*kron(eye(Nz),toeplitz(z,z));

% Solve the system
u = (D2x + D2y)\(f-ub);
% Convert u from a column vector to a matrix to make it easier to work with
% for plotting.
u = reshape(u,Nx,Nz);
% Append on to u the boundary values from the Dirichlet condition.
% u = [[feval(bfunc,x(1,1:n+2),y(1,1:n+2))];...
% [[feval(bfunc,x(2:n+1,1),y(2:n+1,1))] u ...
% [feval(bfunc,x(2:n+1,n+2),y(2:n+1,n+2))]];...
% [feval(bfunc,x(n+2,1:n+2),y(n+2,1:n+2))]];

n2=u;
