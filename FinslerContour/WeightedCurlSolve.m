% Copyright Jean-Marie Mirebeau, 2015. jm(dot)mirebeau(at)gmail(dot)com

% Solves :
% min integral |u1|^2 rho_1, subject to curl u1 =  f
% min integral |u2|^2 rho_2, subject to curl u2 = -f

% Outputs : u1, u2, p such that u1+u2= rot p

% Internally uses WeightedCurlSolveRaw, followed with some averaging to 
% go back to standard grid.

function [u1,u2,p] = WeightedCurlSolve(f,rho1,rho2,gridScale)
    assert(ndims(f)==2);
    assert(all(size(f)==size(rho1)));
    assert(all(size(f)==size(rho2)));

    [u1h,u1v]=WeightedCurlSolveRaw( f,rho1);
    [u2h,u2v]=WeightedCurlSolveRaw(-f,rho2);
    
    p=AveragedPotential(u1h+u2h,u1v+u2v) *gridScale*gridScale;

    u1=AverageOnGrid(u1h,u1v) * gridScale;
    u2=AverageOnGrid(u2h,u2v) * gridScale;
end

function u = AverageOnGrid(uh,uv) 
    n=size(uh,2); m=size(uv,1);
    assert(size(uh,1)==m+1);
    assert(size(uv,2)==n+1);
    
    u=zeros(2,m,n);
    u(1,:,:) = (uh(1:m,:)+uh(2:(m+1),:))/2.;
    u(2,:,:) = (uv(:,1:n)+uv(:,2:(n+1)))/2.;
end

function p = AveragedPotential(uh,uv)
    n=size(uh,2); m=size(uv,1);
    assert(size(uh,1)==m+1);
    assert(size(uv,2)==n+1);
    
    p = [0;cumsum(uv(:,1))] * ones(1,n+1) + [zeros(m+1,1),cumsum(uh,2)];
    p=(p(1:m,1:n)+p(2:(m+1),1:n)+p(1:m,2:(n+1))+p(2:(m+1),2:(n+1)))/4;
end