% Copyright Jean-Marie Mirebeau, 2015. jm(dot)mirebeau(at)gmail(dot)com

% Chan-Vese type segmentation using Finsler fast marching.

% Improvements to do, for faster convergence when f takes large values: 
% - what weight for the curl solve ? 1/d(x,C) does not work in full space. Square of it ?
% - How to modify the metric when omega is large ? Id+omega x omega ? (Nicely consistent in neighborhood.)
%   Or alpha*Id+omega x omega * F(omega^2/alpha), with F(1/2 - .)=0, F(1 + .)=1.

% Optimum should be a circle of radius 2+Sqrt(3) = 3.73... for rho=-4/(1+r2); 
% We find a slightly larger radius unfortunately. Unclear why.


addpath('/Users/mirebeau/bin/2015/Matlab/MFiles') % -> AnisotropicFastMarching_EarlyAbort.mexmaci64
addpath('../GetOptions')

eps=0.05;
[x,y]=meshgrid(-6:eps:6,-7:eps:7);
r2=x.*x+y.*y;
rho=-4./(1+r2);
sigma=ones(size(rho));

options.PathWidth=0.5;
options.x = x;
options.y = y;

a=3; n=ceil(a/eps);
%First path counterclockwise, second path clockwise
path1=[linspace(-a,a,2*n);linspace(0,-a,n),linspace(-a,0,n)];
path2=[linspace(-a,a,2*n);linspace(0,a,n),linspace(a,0,n)];

options.GridScale = eps;
options.Origin = [x(1,1);y(1,1)];

if true
    options.ShowDistanceFromPath=false;
    options.ShowDensities=false;
    options.ShowVectorFields=false;
    options.ShowVectorFieldsNorms=true;
    options.ShowDistanceMaps=false;
    options.ShowImprovedPaths=true;

    pause on;
    for i=1:10
        [path1,path2] = ImprovePaths_Loc(path1,path2,rho,sigma,options);
        %Note: to use the faster, better version ImprovePaths_Loc, a small modification is needed in AnisotropicFastMarching_EarlyAbort
    end;
end
