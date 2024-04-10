% Copyright Jean-Marie Mirebeau, 2015. jm(dot)mirebeau(at)gmail(dot)com

% Solves :
% min integral |u|^2 rho, subject to curl u = f

% If f : mxn grid, then uh : mx(n+1), and uv : (m+1)xn 

function [uh,uv] = WeightedCurlSolveRaw(f,rho)
    assert(ndims(f)==2);
    assert(all(size(f)==size(rho)));    
    m=size(f,1); n=size(f,2);
    
    column = zeros( 1, 8*m*n +m*(n+1) +(m+1)*n );
    row = column;
    coeff = column;
    
    %Numbering of constraints and unknowns
    cstr = 1:(m*n);
    unkh=1:((m+1)*n); unkh=cstr(end)+unkh;
    unkv=1:(m*(n+1)); unkv=unkh(end)+unkv;
    
    %Creating the constraint matrix
    unkh=reshape(unkh,m+1,n);
    unkv=reshape(unkv,m,n+1);
    
    unkh1=reshape(unkh(1:m,    :),[1,m*n]);
    unkh2=reshape(unkh(2:(m+1),:),[1,m*n]);
    unkv1=reshape(unkv(:,    1:n),[1,m*n]);
    unkv2=reshape(unkv(:,2:(n+1)),[1,m*n]);
    
    pos=1:(4*m*n);
    unit = ones([1,m*n]);
    column(pos) = [unkh1,unkh2,unkv1,unkv2];
    row(pos)    = [cstr,cstr,cstr,cstr];
    coeff(pos)  = [unit,-unit,-unit,unit];
        
    % Creating the Lagrange multiplier matrix, transpose of constraint
    pos2=pos(end)+pos;
    column(pos2)=row(pos);
    row(pos2)=column(pos);
    coeff(pos2)=coeff(pos);
    
    %Weights
    unkh=unkh(:)'; 
    unkv=unkv(:)';
    pos3=(pos2(end)+1):size(column,2);
    column(pos3) = [unkh,unkv];
    row(pos3) = [unkh,unkv];
    rhoh=imfilter(rho,[1;1]/2.,'replicate','full');
    rhov=imfilter(rho,[1,1]/2.,'replicate','full');
    coeff(pos3)=[rhoh(:)',rhov(:)'];

    % Rhs and matrix solve
    
    unkn = m*n+m*(n+1)+(m+1)*n;
    rhs = [f(:);zeros(m*(n+1)+(m+1)*n,1)];
    A=sparse(column, row, coeff, unkn, unkn);
    %opts.SYM = true;
    %sol = linsolve(A,rhs,opts);
    sol = A\rhs;
    
    % Output
    uh=sol(unkh);
    uv=sol(unkv);
    uh=reshape(uh,m+1,n);
    uv=reshape(uv,m,n+1);
    
end