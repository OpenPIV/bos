function [n2] = Poisson_directMod(Nx, Nz, dxb, dxc, dzb, dzc,rhs)
idx=reshape(1:(Nx+2)*(Nz+2), Nx+2, Nz+2);

ii=zeros(5*Nx*Nz,1);
jj=zeros(5*Nx*Nz,1);
vv=zeros(5*Nx*Nz,1);

% Compute the terms of the matrix [L] using a stencil of 5 points
a=1; ss=5; 
for k=2:Nz+1
  for i=2:Nx+1
    ii(a:a+ss-1) = idx(i, k) * ones(1,ss);
    jj(a:a+ss-1) = [ idx(i-1,k), idx(i,k-1), idx(i,k), ...
                     idx(i,k+1), idx(i+1,k)];
    vv(a:a+ss-1) = [ 1/(dxc(i)*dxb(i-1)) ...
                   , 1/(dzc(k)*dzb(k-1)) ...
                   , - 1/(dxc(i)*dxb(i-1)) - 1/(dxc(i)*dxb(i)) ...
                     - 1/(dzc(k)*dzb(k-1)) - 1/(dzc(k)*dzb(k)) ...
                   , 1/(dzc(k)*dzb(k)) ...
                   , 1/(dxc(i)*dxb(i))];
    a            = a + ss;
  end
end
L = sparse(ii, jj, vv, (Nx+2)*(Nz+2),(Nx+2)*(Nz+2));

a=1;
ss = 2;
for i=1:Nx+2
  ii(a:a+1) = idx(i, 1) * ones(1,ss); 
  jj(a:a+1) = [idx(i,1) idx(i,2)]; 
  vv(a:a+1) = [-1 1];
  a = a + 2;
end

for i=1:Nx+2
  ii(a:a+1) = idx(i, end) * ones(1,ss); 
  jj(a:a+1) = [idx(i,end-1) idx(i,end)]; 
  vv(a:a+1) = [-1 1];
  a = a + 2;
end

for k=2:Nz+1
  ii(a:a+1) = idx(1, k) * ones(1,ss); 
  jj(a:a+1) = [idx(1, k) idx(2,k)]; 
  vv(a:a+1) = [-1 1];
  a = a + 2;
end

for k=2:Nz+1
  ii(a:a+1) = idx(end, k) * ones(1,ss); 
  jj(a:a+1) = [idx(end-1,k) idx(end,k)]; 
  vv(a:a+1) = [-1 1];
  a = a + 2;
end

Bp = sparse(ii(1:a-1), jj(1:a-1), vv(1:a-1),(Nx+2)*(Nz+2),(Nx+2)*(Nz+2));

% Boundary conditions
%Bp(idx(1:end,1), idx(1,1))=1.43; 
% Bp(idx(1:end,end), idx(1,end))=1; 

L = L + Bp;

p2 = reshape(L \ rhs(:), Nx+2, Nz+2);
pg2=p2/(max(max(p2))-min(min(p2)));
PG2_gray=pg2*255;
n_max=1.43;    
n_min=1.332;
n2= scaledata(PG2_gray,n_min,n_max);

