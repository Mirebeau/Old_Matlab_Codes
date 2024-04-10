function data = NonDivergenceOpt(data,action)
% Copyright Jean-Marie Mirebeau, 2016, CNRS, University Paris-Sud
% This function implements the operator
%    inf_i Tr(D_i Hess(u)) - q_i)

% There are four call possibilities, depending on the parameter action
% -- 'Init'
%   Input : data.nonDivOpt ->
%       - tensors. Format : [SymDimension,nPoints,nTensors]
%   Output : data.nonDivOpt ->
%       - reqOffsets2. Format : [Dimension,nReqOffsets2]. Symmetric finite differences needed by the scheme.
%   Internal : data.nonDIvOpt ->
%       - creates offsets2 [Dimension,nPoints,SymDimension,nTensor], weights2 [1,nPoints,SymDimension,nTensor]. 
%           Tensors decomposition based on obtuse superbases (for internal use).

% -- 'OffsetIndex'
%   Input : data ->
%       - offsets2. Format [Dimension,nOffsets2]. Symmetric finite differences provided.
%   Internal : data.nonDivOpt ->
%       - creates offsets2Index [1,nPoints,SymDimension,nTensor]
%       - deletes reqOffsets2,  offsets2

% -- 'Value'
%   Input :  
%       - data.diff2 [1,nPoints,nOffsets2]. Second order finite differences.
%       - data.nonDivOpt.constantTerm [1,nPoints].
%       - data.nonDivOpt.maximize (bool, default false). Operator sup instead of inf. 
%   Output : data.nonDivOpt ->
%       - values [1,nPoints]
%       - optIndex [1,nPoints]. Active tensor field.
%       - validGradient (bool). Set to false.

% -- 'Gradient'
%   Output : 
%       - grad2 sparse(nPoints,nRefOffsets2*nPoints).
%       - grad1 sparse(nPoints,nRefOffsets1*nPoints). Null for this operator.
%       - grad0 sparse(nPoints,nPoints). Null for this operator.
%       - validGradient (bool). Set to true.
%

assert(isfield(data,'nonDivOpt'));

if strcmp(action,'Init')
    data = NonDivergenceOpt_Init(data);
elseif strcmp(action,'OffsetIndex') 
    data = NonDivergenceOpt_OffsetIndex(data);
elseif strcmp(action,'Value')
    data = NonDivergenceOpt_Value(data);
elseif strcmp(action,'Gradient')
    data = NonDivergenceOpt_Gradient(data);
else 
    disp('NonDivergenceOpt : unrecognized action')
end

end

% ---------------------------------------

function data = NonDivergenceOpt_Init(data)
tensorsSize = size(data.nonDivOpt.tensors);  
SymDimension = tensorsSize(1); % Should we call this stencil dimension ?
nPoints = tensorsSize(2);
nTensors = tensorsSize(3);
Dimension = floor(sqrt(2*SymDimension));

assert(all(tensorsSize==[SymDimension,nPoints,nTensors]));
offsets = zeros([Dimension,nPoints,SymDimension,nTensors]);
weights = zeros([1,nPoints,SymDimension,nTensors]);

for iTensor = 1:nTensors
    [offsets(:,:,:,iTensor),weights(:,:,:,iTensor)] = Sym_Decomposition(data.nonDivOpt.tensors(:,:,iTensor));
end

%TODO : normalize offsets so that first non-zero component is positive (reduces a bit memory usage)
offsets = reshape(Offset_Normalize(reshape(offsets,[Dimension,nPoints*SymDimension*nTensors])),[Dimension,nPoints,SymDimension,nTensors]);

data.nonDivOpt.offsets2 = offsets;
data.nonDivOpt.weights2 = weights;

offsets = reshape(offsets,[Dimension,nPoints*SymDimension*nTensors]);
data.nonDivOpt.reqOffsets2 = union(offsets',offsets(:,1)','rows')';
end

% ---------------------------------------

function data = NonDivergenceOpt_OffsetIndex(data)

offsets2 = data.nonDivOpt.offsets2; %size is [Dimension,nPoints,SymDimension,nTensor]
Dimension = size(offsets2,1);
nPoints = size(offsets2,2);
SymDimension = size(offsets2,3);
nTensors = size(offsets2,4);

refOffsets2 = data.offsets2;
nRefOffsets2 = size(refOffsets2,2);

refOffsets2 = reshape(refOffsets2,[Dimension,nRefOffsets2]);
offsets2=reshape(offsets2,[Dimension,nPoints*SymDimension*nTensors]);
[found,offsets2Index] = ismember(offsets2',refOffsets2','rows');
assert(all(found));

% Give absolute index
ptIndex = repmat(1:nPoints,[1,SymDimension*nTensors]);
offsets2Index = sub2ind([nPoints,nRefOffsets2],ptIndex,offsets2Index');

offsets2Index = reshape(offsets2Index,[1,nPoints,SymDimension,nTensors]);

data.nonDivOpt.offsets2Index = offsets2Index;

rmfield(data.nonDivOpt, 'offsets2');
rmfield(data.nonDivOpt, 'reqOffsets2');

end

% ---------------------------------------

function data = NonDivergenceOpt_Value(data)
tensorSize = size(data.nonDivOpt.tensors);
SymDimension = tensorSize(1);
nPoints = tensorSize(2);
nTensors = tensorSize(3);

if isfield(data.nonDivOpt,'maximize') && data.nonDivOpt.maximize
    values = -Inf*ones([1,nPoints]);
else 
    values = Inf*ones([1,nPoints]);    
end

optIndex = zeros([1,nPoints]);
for iTensor = 1:nTensors
    offsets2Index = data.nonDivOpt.offsets2Index(1,:,:,iTensor);
    weights = data.nonDivOpt.weights2(1,:,:,iTensor);
    
    diff2 = reshape(data.diff2(offsets2Index(:)),[1,nPoints,SymDimension]);    
    tValues = sum(diff2.*weights,3) - data.nonDivOpt.constantTerm(1,:,iTensor);
    
    if isfield(data.nonDivOpt,'maximize') && data.nonDivOpt.maximize
        pos = tValues>values;
    else
        pos = tValues<values; % exludes nan
    end
    values(pos) = tValues(pos);
    optIndex(pos) = iTensor;
end

data.nonDivOpt.values = values;
data.nonDivOpt.optIndex = optIndex;
data.nonDivOpt.validGradient = false;

end

% ---------------------------------------

function data = NonDivergenceOpt_Gradient(data)
assert(data.nonDivOpt.validGradient==false);
tensorSize = size(data.nonDivOpt.tensors);
SymDimension = tensorSize(1);
nPoints = tensorSize(2);
nTensors = tensorSize(3);

optIndex = data.nonDivOpt.optIndex;

% Get the appropriate second order difference indices
ptIndex = repmat(1:nPoints,[SymDimension,1]); ptIndex=ptIndex(:);
o2Index = repmat(1:SymDimension,[1,nPoints]); o2Index=o2Index(:);
optIndex = repmat(optIndex,[SymDimension,1]); optIndex=optIndex(:);

ind = sub2ind([nPoints,SymDimension,nTensors],ptIndex,o2Index,optIndex);

Rows = ptIndex;
Cols = data.nonDivOpt.offsets2Index(ind);
Vals = data.nonDivOpt.weights2(ind);

nRefOffsets2 = size(data.offsets2,2);
data.nonDivOpt.grad2 = sparse(Rows,Cols,Vals,nPoints,nRefOffsets2*nPoints);

data.nonDivOpt.grad0 = sparse([nPoints,nPoints]);
nRefOffsets1 = size(data.offsets1,2);
data.nonDivOpt.grad1 = sparse([nPoints,nRefOffsets1*nPoints]);
data.nonDivOpt.validGradient = true;

end