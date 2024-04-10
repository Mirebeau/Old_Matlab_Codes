% Copyright Jean-Marie Mirebeau, 2015. jm(dot)mirebeau(at)gmail(dot)com

% Solves :
% min integral |u|^2 rho, subject to curl u = f

% If f : mxn grid, then uh : mx(n+1), and uv : (m+1)xn 

function [uh,uv] = WeightedCurlSolveRaw_Loc(f,rho)
    assert(ndims(f)==2);
    assert(all(size(f)==size(rho)));
    m=size(f,1); n=size(f,2);

        rhoh=imfilter(rho,[1;1]/2.,'replicate','full');
    rhov=imfilter(rho,[1,1]/2.,'replicate','full');

    
    %Zero entries in rho are for mask
    mask = rho>0;
    cstrLoc = find(mask>0); %Indices in array
    cstr = (1:numel(cstrLoc))';
    
    rhoh=imfilter(mask,[1;1]/2.,'replicate','full'); 
    rhov=imfilter(mask,[1,1]/2.,'replicate','full'); 
    
    uhLoc=find(rhoh>0);
    uvLoc=find(rhov>0);
    
    uhUnk=zeros([m+1,n]);
    uvUnk=zeros([m,n+1]);
    uhUnk(uhLoc) = cstr(end)+(1:numel(uhLoc));
    uvUnk(uvLoc) = uhUnk(uhLoc(end))+(1:numel(uvLoc));
        
    % Creating the constraint matrix
    column=zeros(8*numel(cstrLoc)+numel(uhLoc)+numel(uvLoc),1);
    row=column;
    coeff=column;
    
    pos=1:(4*numel(cstr));
    unit=ones([numel(cstr),1]);
    
    uhUnk1=uhUnk(1:m,       :); uhUnk1=uhUnk1(cstrLoc);
    uhUnk2=uhUnk(2:(m+1),   :); uhUnk2=uhUnk2(cstrLoc);
    uvUnk1=uvUnk(:,       1:n); uvUnk1=uvUnk1(cstrLoc);
    uvUnk2=uvUnk(:,   2:(n+1)); uvUnk2=uvUnk2(cstrLoc);
    
    column(pos) = [uhUnk1;uhUnk2;uvUnk1;uvUnk2];
    row(pos)=[cstr;cstr;cstr;cstr];
    coeff(pos)=[unit;-unit;-unit;unit];
    

    % Creating the Lagrange multiplier matrix, transpose of constraint
    pos2=pos(end)+pos;
    column(pos2)=row(pos);
    row(pos2)=column(pos);
    coeff(pos2)=coeff(pos);
    
    %Weights
    pos3=(pos2(end)+1):numel(column);
    'WeightedCurlSolveRaw, sizes...'
    size(uhUnk(uhLoc))
    size(pos3)
    size(cstrLoc)
    size(uhLoc)
    size(uvLoc)
    column(pos3) = [uhUnk(uhLoc);uvUnk(uvLoc)];
    row(pos3) = column(pos3);
    coeff(pos3)=[rhoh(uhLoc);rhov(uvLoc)];

    % Rhs and matrix solve
    
    unkn = numel(cstr)+numel(uhLoc)+numel(uvLoc);
    rhs = [f(cstrLoc);zeros(numel(uhLoc)+numel(uvLoc),1)];
    A=sparse(column, row, coeff, unkn, unkn);
    %opts.SYM = true;
    %sol = linsolve(A,rhs,opts);
    sol = A\rhs;
    
    % Output
    uh=zeros([m+1,n]);
    uv=zeros([m,n+1]);
    uh(uhLoc)=sol(numel(cstrLoc)+(1:numel(uhLoc)));
    uv(uvLoc)=sol(numel(cstrLoc)+numel(uhLoc)+(1:numel(uvLoc)));    
end