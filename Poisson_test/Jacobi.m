function [n2]=Jacobi(Nx, Nz,Rhs)


tol=1e-6;   % tollerance
err=1;      % error
k=0;        % iteration counter

% Initialize the n-vale
nki=Rhs;

% Output column heading
 %fprintf(' k |')
 for j=1:Nz-2
        for i=1:Nx-2 
     %fprintf(' n(%1i,%1i) |', j,i)
        end
 end
 
     %fprintf(' error\n')

% Iterative Jacobi untill convergence
while err>tol
    
   k=k+1;
% loop through computational node (inside the matrix)

    for j=2:Nx
     for i=2:Nz
         nki(j,i)=0.25*(Rhs(j+1,i)+Rhs(j,i+1)+Rhs(j-1,i)+Rhs(j,i-1));
    % fprintf('%8.6f|',nki(j,i))
     end
    end

    % calculate the error
    
    err=sqrt(sum(sum(nki-Rhs).^2));
  %  fprintf('%8.6f|',err)

    % update n
    Rhs=nki;
end
n2=nki;
