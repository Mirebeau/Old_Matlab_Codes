% Method parameters
Dimension=2;
n=99;
gridScale = 1./n;
eps=4*gridScale;
nIter = 30;
seed = [0.5,0.5];
xCoords = gridScale*((1:n)-0.5);
yCoords = gridScale*((1:(n+1))-0.5);

% Finsler metric

dims = [numel(xCoords),numel(yCoords)];
nPoints=prod(dims);
SymDimension=Dimension*(Dimension+1)/2;
m=ones([SymDimension,nPoints]);
m(2,:)=0.;
w=0.*ones([Dimension,nPoints]);
w(1,:)=0.2;
w(2,:)=0.3;

% Import required libs

addpath('../NumericalScheme');
addpath('../../ParallelMatrix');

% Define differential operator

[s,v,a]=EikoFinslerLinear(m,w,eps);

% Build its sparse matrix

options.dims=dims;
options.gridScale=gridScale;
[Rows,Cols,Vals] = DifferentialOperatorMatrix(s,v,a,options);

A=sparse(Rows,Cols,Vals,nPoints,nPoints);

% Solve to obtain exp(-distance/eps)
b=zeros(dims);
b(ceil(n*seed(1)),ceil(n*seed(2)))=1;
b=reshape(b,[nPoints,1]);
sol = A\b;
approxDist = -eps*log(sol);
approxDist = approxDist - min(approxDist(:));
approxDist=reshape(approxDist,dims);
imagesc(approxDist);
pause;

[x,y] = meshgrid(xCoords,yCoords); x=x'; y=y';
u = [x(:)-seed(1),y(:)-seed(2)]; u=u';
exactDist = SymVecVec_EvalFinsler(m,w,u); 
exactDist=reshape(exactDist,dims);
imagesc(exactDist);
pause;

imagesc(exactDist-approxDist);
pause;

% Show a one dimensional slice through the seed
plot(xCoords,exactDist(:,ceil(n*seed(2))),xCoords,approxDist(:,ceil(n*seed(2))));

% Compute the geodesic directions, which should be straight to origin in this case.
options.eps=eps;
sol=reshape(sol,[1,nPoints]);
dir=EikoFinslerDirection(s,v,sol,options);

dx=reshape(dir(1,:),dims); dy=reshape(dir(2,:),dims);
quiver(x',y',dx',dy')