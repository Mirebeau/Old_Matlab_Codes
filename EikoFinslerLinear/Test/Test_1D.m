% Method parameters
Dimension=1;
n=199;
gridScale = 1./n;
eps=4*gridScale;
nIter = 30;
seed = 0.6;

% Finsler metric

nPoints=n^Dimension;
SymDimension=Dimension*(Dimension+1)/2;
m=ones([SymDimension,nPoints]);
w=-0.1*ones([Dimension,nPoints]);

% Import required libs

addpath('../NumericalScheme');
addpath('../../ParallelMatrix');

% Define differential operator

[s,v,a]=EikoFinslerLinear(m,w,eps);

% Build its sparse matrix

options.dims=n;
options.gridScale=gridScale;
[Rows,Cols,Vals] = DifferentialOperatorMatrix(s,v,a,options);

A=sparse(Rows,Cols,Vals,nPoints,nPoints);

% Solve to obtain exp(-distance/eps)
b=zeros([nPoints,1]);
b(ceil(n*seed))=1;
sol = A\b;

xAxis=((1:n)-0.5)*gridScale;
yAxis=-eps*log(sol); yAxis=yAxis-min(yAxis(:));
shiftAxis = xAxis-seed;
plot(xAxis,yAxis,xAxis,SymVecVec_EvalFinsler(m,w,shiftAxis)); %abs(shiftAxis)*m(1)+shiftAxis*w(1)
pause;

% Iterate to obtain a transported heat kernel
for iter=2:nIter
    sol = A\sol;
end

diffusionTime = nIter*eps^2;
transportTime = 2*nIter*eps;

xAxis=(1:n)*gridScale;
yAxis=-diffusionTime*log(sol); yAxis=yAxis-min(yAxis(:));
shiftAxis = xAxis-(seed+w(1)*transportTime);
plot(xAxis,yAxis,xAxis,(1/4)*shiftAxis.^2*m(1));
