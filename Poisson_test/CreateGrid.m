function [xc zc dxb dxc dzc dzb Nx Nz]= CreateGrid(Lx, Lz,Rx,Ry)
% 
Nx=Rx-2;
Nz=Ry-2;

xi=linspace(0,Lx,Nx+1);
xb = xi;
dxc=ones(length(xi));
dxb=ones(length(xi));

xc=linspace(0,Lx,Nx+2);
zi=linspace(0,Lz,Nz+1);
zb = zi;
zc=linspace(0,Lz,Nz+2);
dzc=ones(length(zi));
dzb=ones(length(zi));



