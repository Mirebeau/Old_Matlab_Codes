% A 2D test, with a variable finsler metric, and an exact radial solution
% Based on Proposition 4.2 of J.-M. Mirebeau, Efficient Fast Marching with Finsler Metrics, Numerische Math

% Method parameters
Dimension=2;
n=199;
gridScale = 1./n;
eps=4*gridScale;
nIter = 30;
seed = [0.5,0.5];

% Exact solution

%radialDistDeriv = @(r) 0.5*(1+1./(1+r));
%radialDist = @(r) 0.5*(r+log(1+r));

radialDistDeriv = @(r) min(1,0.5*(1+(4*r-1).^2));
radialDist = @(r) 0.5*(r+(4*r-1).^3/12);

coords = gridScale*((1:n)-0.5);
[x,y] = meshgrid(coords,coords); x=x'; y=y';
r = sqrt((x-seed(1)).^2.+(y-seed(2)).^2);
exactDist  = radialDist(r);
imagesc(exactDist);
pause;

% Finsler metric

radialDensity = @(r) sqrt(1-radialDistDeriv(r).^2);
nPoints=n^Dimension;
SymDimension=Dimension*(Dimension+1)/2;
m=ones([SymDimension,nPoints]);
m(2,:)=0.;
w=zeros([Dimension,nPoints]);

w(1,:) = reshape(-radialDensity(r).*(y-seed(2))./r,[1,nPoints]);
w(2,:) = reshape( radialDensity(r).*(x-seed(1))./r,[1,nPoints]);

w(:,r==0)=0;

% Import required libs

addpath('../NumericalScheme');
addpath('../../ParallelMatrix');

% Define differential operator

[s,v,a]=EikoFinslerLinear(m,w,eps);

% Build its sparse matrix

dims=[n,n];
options.dims=dims;
options.gridScale=gridScale;
[Rows,Cols,Vals] = DifferentialOperatorMatrix(s,v,a,options);

A=sparse(Rows,Cols,Vals,nPoints,nPoints);

% Solve to obtain exp(-distance/eps)
b=zeros(dims);
b(ceil(n*seed(1)),ceil(n*seed(2)))=1;
b=reshape(b,[nPoints,1]);
sol = A\b;
sol=reshape(sol,[n,n]);
approxDist = -eps*log(sol);
approxDist = approxDist - min(approxDist(:));
imagesc(approxDist);
pause;

imagesc(exactDist-approxDist);
pause;

% Show a one dimensional slice through the seed
ySlice = ceil(n*seed(2));
exactSlice = exactDist(:,ySlice);
approxSlice = approxDist(:,ySlice);
plot(coords,exactSlice,coords,approxSlice);

plot(coords(2:n),diff(exactSlice)/gridScale,coords(2:n),diff(approxSlice)/gridScale);
pause;

% Show gradient

options.eps = eps;
sol=reshape(sol,[1,nPoints]);
dir=EikoFinslerDirection(s,v,sol,options);
dx=reshape(dir(1,:),dims); dy=reshape(dir(2,:),dims);
quiver(x',y',dx',dy')
