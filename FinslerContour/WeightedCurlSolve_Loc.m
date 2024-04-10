% Copyright Jean-Marie Mirebeau, 2015. jm(dot)mirebeau(at)gmail(dot)com

% Solves :
% min integral |u1|^2 rho_1, subject to curl u1 =  f
% min integral |u2|^2 rho_2, subject to curl u2 = -f

% Outputs : u1, u2, p such that u1+u2= rot p

% Internally uses WeightedCurlSolveRaw, followed with some averaging to 
% go back to standard grid.

function [u1,u2,p] = WeightedCurlSolve_Loc(f,rho1,rho2,gridScale,rhop)
    assert(ndims(f)==2);
    assert(all(size(f)==size(rho1)));
    assert(all(size(f)==size(rho2)));

    [u1h,u1v]=WeightedCurlSolveRaw_Loc( f,rho1);
    [u2h,u2v]=WeightedCurlSolveRaw_Loc(-f,rho2);
    
    p=AveragedPotential(u1h+u2h,u1v+u2v,rhop) *gridScale*gridScale;

    u1=AverageOnGrid(u1h,u1v) * gridScale;
    u2=AverageOnGrid(u2h,u2v) * gridScale;
end

function u = AverageOnGrid(uh,uv) 
    n=size(uh,2); m=size(uv,1);
    assert(size(uh,1)==m+1);
    assert(size(uv,2)==n+1);
    
    u=zeros(2,m,n); %Sparse tensors not supported...
    u(1,:,:) = (uh(1:m,:)+uh(2:(m+1),:))/2.;
    u(2,:,:) = (uv(:,1:n)+uv(:,2:(n+1)))/2.;
end

% Problem here ??
function p = AveragedPotential(uh,uv,rhop)
    m=size(rhop,1); n=size(rhop,2);
    assert(all(size(uh)==[m+1,n]));
    assert(all(size(uv)==[m,n+1]));
    
    % Unknowns locations and numbering
    mask=imfilter(rhop,[1,1;1,1]/4.,'replicate','full');
    pLoc = find(mask>0);
    pUnk = zeros([m+1,n+1]);
    pUnk(pLoc) = 1:numel(pLoc);
    
    hMask = imfilter(rhop,[1;1]/2,'replicate','full');
    vMask = imfilter(rhop,[1,1]/2,'replicate','full');
    hLoc = find(hMask>0);
    vLoc = find(vMask>0);
    
    pUnkH1 = pUnk(1:(m+1),     1:n);
    pUnkH2 = pUnk(1:(m+1), 2:(n+1));
    pUnkV1 = pUnk(    1:m, 1:(n+1));
    pUnkV2 = pUnk(2:(m+1), 1:(n+1));
    
    % Creating the sparse matrix and solving
    
    hN = numel(hLoc); vN=numel(vLoc);
    column = [1:hN,        1:hN,        hN+(1:vN),   hN+(1:vN),    hN+vN+1];
    row =    [pUnkH1(hLoc);pUnkH2(hLoc);pUnkV1(vLoc);pUnkV2(vLoc); 1];
    coeff =  [-ones(1,hN), ones(1,hN),  -ones(1,vN), ones(1,vN),   1];
    
    rhs = [uh(hLoc);uv(vLoc); 0];
    A=sparse(column,row,coeff);
    'AveragedPotential, sizes...'
    size(rhs)
    size(A)
    
    sol = A\rhs;
    p = zeros([m+1,n+1]);
    p(pLoc)=sol;
    
    %Averaging
    p=(p(1:m,1:n)+p(2:(m+1),1:n)+p(1:m,2:(n+1))+p(2:(m+1),2:(n+1)))/4;
end