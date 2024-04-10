addpath('../NumericalScheme');
addpath('../../ParallelMatrix');

clear('data')
dims=[9,10];
nPoints=prod(dims);
Dimension = size(dims,2);
SymDimension = Dimension*(Dimension+1)/2;

data.nonDivOpt.tensors = ones([SymDimension,nPoints,2]);
data.nonDivOpt.tensors(2,:,1)=0;
data.nonDivOpt.tensors(2,:,2)=0.5;

nTensors = size(data.nonDivOpt.tensors,3);
data=NonDivergenceOpt(data,'Init'); % Diffusion tensors decomposition

data.geometry.offsets = union(data.nonDivOpt.reqOffsets2',-data.nonDivOpt.reqOffsets2','rows')';
data.geometry.uniformOffsets=true;
data=GridCrop(dims,data); % Construction of the box indices

data.offsets1=zeros([Dimension,0]);
data.offsets2=data.nonDivOpt.reqOffsets2;
data.boundaryValue = @(x) x(1,:).^2+x(2,:).^2;
data=DifferenceMatrix(data);

data=NonDivergenceOpt(data,'OffsetIndex'); 

u = data.boundaryValue(data.geometry.points)';
data.diff0 = u; % zeroth order differences are alues
data.diff1 = data.linear.A1*u+data.linear.b1; 
data.diff2 = data.linear.A2*u+data.linear.b2;

nTensors = size(data.nonDivOpt.tensors,3);
data.nonDivOpt.constantTerm = zeros([1,nPoints,nTensors]);

data=NonDivergenceOpt(data,'Value'); % Computation of the values
data=NonDivergenceOpt(data,'Gradient'); % Computation of the gradient