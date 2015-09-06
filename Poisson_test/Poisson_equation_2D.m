clc
clear all
close all

% Solving the 2-D Poisson equation by the Finite Difference
...Method 
% Numerical scheme used is a second order central difference in space
...(5-point difference)
    
[Mconversion,Const,Lx,Lz,val_up,val_down,nx_pixel,ny_pixel,overlap_x,overlap_y]=Parameters();

% Load the displacement only for debug

%load Displacement_POisson.mat
load Displacement_POisson1.mat
% 
Displ=Displacement_POisson;
[Displ]=crop_field(Displacement_POisson,Lx,Lz);
[Rhs]=create_RHS(Displ);

[Nx,Nz]=size(Rhs);

%%
%Specifying parameters (check why It does not work if we change them)

niter=100000;                      %Number of iterations 
dx=Lx/(Nx-1);                     %Width of space step(x)
dy=Lz/(Nz-1);                     %Width of space step(y)
x=0:dx:Lx;                        %Range of x(0,2) and specifying the grid points
y=0:dy:Lz;                        %Range of y(0,2) and specifying the grid points
b=zeros(Nx,Nz);                  %Preallocating b
pn=zeros(Nx,Nz);                 %Preallocating pn

%%
% Initial Conditions
p=zeros(Nx,Nz);                  %Preallocating p

%% 
%Rhs=Const.*fliplr(Rhs);   

%b=(Rhs)/max(max(Rhs));
%b=(Rhs')/norm(Rhs');
b=(Rhs);


%%
i=2:Nx-1;
j=2:Nz-1;

% Poisson equation solution (iterative) method

tol = 1e-4;			 % error is 1%
maxerr = inf;	 % initial error
iter = 0;
pn=p;

while maxerr > tol
    iter = iter + 1;
    disp(['Iteration no. ',num2str(iter)]);

    
%Explicit iterative scheme with C.D in space (5-point difference)
%for it=1:niter
   
    p(i,j)=((dy^2*(pn(i+1,j)+pn(i-1,j)))+(dx^2*(pn(i,j+1)+pn(i,j-1)))-(b(i,j)*dx^2*dy*2))/(2*(dx^2+dy^2));
    %Boundary conditions 
    
   % Neumann's conditions   % dp/dx|end=dp/dx|end-1 
    p(1,:)=p(2,:);
    p(end,:)=p(end-1,:);
    
    
%   %  Dirichlet's conditions 
%     p(:,1)=1.43;
%     p(:,end)=1.32;
    
   % Neumann's conditions
    p(:,1)=p(:,2);
    p(:,end)=p(:,end-1);

    maxerr = max(max(abs((p-pn)./p)));
    disp(['Maximum error is  ',num2str(maxerr)]);
    pn=p;
   
end			% as long the error larger than tolerance, continue
%%
%Plotting the solution

PG2_gray=p*255;
n_max=1.43;    
n_min=1.332;
n2= scaledata(PG2_gray,n_min,n_max);

figure 
subplot(121)
h=surf(x,y,p','EdgeColor','none');       
shading interp
%axis([-0.5 2.5 -0.5 2.5 -100 100])
%title({'2-D Poisson equation';['{\itNumber of iterations} = ',num2str(it)]})
xlabel('Spatial co-ordinate (x) \rightarrow')
ylabel('{\leftarrow} Spatial co-ordinate (y)')
zlabel('Solution profile (P) \rightarrow')

subplot(122)
h=surf(x,y,n2','EdgeColor','none');       
shading interp
%axis([-0.5 2.5 -0.5 2.5 -100 100])
%title({'2-D Poisson equation';['{\itNumber of iterations} = ',num2str(it)]})
xlabel('Spatial co-ordinate (x) \rightarrow')
ylabel('{\leftarrow} Spatial co-ordinate (y)')
zlabel('Solution profile (P) \rightarrow')


% Plot a scaled profile 
figure
plot(n2(5,:),y)