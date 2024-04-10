addpath('../NumericalScheme');
addpath('../../ParallelMatrix');

clear('data')
dims=9;
nPoints=prod(dims);
data.geometry.offsets = [1,-1,2,-2];
data = GridCrop(dims,data);

data.boundaryValue = @(x) (x-5).^2;
data.offsets1 = [1,-2];
data.offsets2 = [1];
data = DifferenceMatrix(data);


data.nonDivOpt.tensors = reshape([ones([1,dims(1)]), 2.*(1:nPoints)/(nPoints+1)],[1,nPoints,2]);
data=NonDivergenceOpt(data,'Init'); % Diffusion tensors decomposition
data=NonDivergenceOpt(data,'OffsetIndex'); % Identification of required offsets

u = data.boundaryValue(data.geometry.points)';
data.diff0 = u; % zeroth order differences are values
data.diff1 = data.linear.A1*u+data.linear.b1; 
data.diff2 = data.linear.A2*u+data.linear.b2;

nTensors = size(data.nonDivOpt.tensors,3);
data.nonDivOpt.constantTerm = zeros([1,nPoints,nTensors]);

data=NonDivergenceOpt(data,'Value'); % Computation of the values
data=NonDivergenceOpt(data,'Gradient'); % Computation of the gradient
