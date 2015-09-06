function [Rhs,Nx,Nz]=create_RHS(Displ)

% This function creates the Rhs of the Poisson's equation. 
% The Rhs is obtained summing dx+dy (defined as Displ.u and Displ.v).

u=Displ.u;
w=Displ.v;
x=Displ.x;
z=Displ.y;

[width,height]=size(u);

% Preallocation for the variables du and dw;
du=zeros(width-2,height-2);
dw=zeros(width-2,height-2);

% Compute the central differences 
for k=2:width-1
for j=2:height-1
    
    du(k,j)=(u(k+1,j)-u(k-1,j))/2*(abs(x(2,2)-x(2,1)));
    dw(k,j)=(w(k,j+1)-w(k,j-1))/2*(abs(z(1,1)-z(1,2)));

end
end

Rhs=(du+dw);
Rhs(isnan(Rhs)) = 0 ;

[Nz,Nx]=size(Rhs);

% create_RHS can be replaced using simply the function gradient()... 