function data = TwoHessian(data,action)
% Copyright Jean-Marie Mirebeau, 2016, CNRS, University Paris-Sud

% This function implements the operator
%    inf_i { Tr(D Hess(u)) - rho sqrt(q(D)) ); Tr(D)=1, D is sb_i acute}
% where (sb_i)_{1<=i <=n} is a collection of superbases,
% and q(D) = Tr(D)^2/(d-1) - Tr(D^2)

% In dimension 2, this amounts to solving a Monge-Ampere equation.
% The coefficient rho is assumed to be non-negative.

% (The discretization is 2D only. (?))

% There are four call possibilities, depending on the parameter action
% -- 'Init'
%   Input : data.twoHessian ->
%       - sb. Format : [Dimension,superbasesN,Dimension+1] the superbases to be used
%   Output : data.twoHessian ->
%       - reqOffsets2. Format : [Dimension,nReqOffsets2]. Symmetric finite differences needed by the scheme.

% -- 'OffsetIndex'
%   Input : data ->
%       - offsets2. Format [Dimension,nOffsets2]. Symmetric finite differences provided.
%   Internal : data.twoHessian ->
%       - creates offsets2Index [1,1,SymDimension,nSuperbases]
%       - deletes reqOffsets2,  offsets2

% -- 'Value'
%   Input :  
%       - data.diff2 [1,nPoints,nOffsets2]. Second order finite differences.
%       - data.twoHessian.rho [1,nPoints].
%   Output : data.nonDivOpt ->
%       - values [1,nPoints]
%       - optIndex [2,nPoints]. Active inf.
%       - validGradient (bool). Set to false.

% -- 'Gradient'
%   Output : 
%       - grad2 sparse(nPoints,nRefOffsets2*nPoints).
%       - grad1 sparse(nPoints,nRefOffsets1*nPoints). Null for this operator.
%       - grad0 sparse(nPoints,nPoints). Null for this operator.
%       - validGradient (bool). Set to true.
%

assert(isfield(data,'twoHessian'));

if strcmp(action,'Init')
    data = TwoHessian_Init(data);
elseif strcmp(action,'OffsetIndex') 
    data = TwoHessian_OffsetIndex(data);
elseif strcmp(action,'Value')
    data = TwoHessian_Value(data);
elseif strcmp(action,'Gradient')
    data = TwoHessian_Gradient(data);
else 
    disp('TwoHessian : unrecognized action')
end

end

% ---------------------------------------

function data = TwoHessian_Init(data)
sb = data.twoHessian.sb; %Superbases to be used
sbSize = size(sb);
Dimension = sbSize(1);
sbN = sbSize(2);
assert(all(sbSize==[Dimension,sbN,Dimension+1]));
SymDimension = Dimension*(Dimension+1)/2;

sbV = zeros([Dimension,sbN,SymDimension]); % offsets associated to the superbases
k=1;
for i=1:Dimension
    for j=(i+1):Dimension
        sbV(:,:,k) = Vec_MultilinPerp(sb(:,:,[1:(i-1),(i+1):(j-1),(j+1):Dimension]));
        k=k+1;
    end
end
assert(k==SymDimension+1);

data.twoHessian.sbV=sbV;
offsets=reshape(sbV,[Dimension,sbN*SymDimension]);
data.twoHessian.reqOffsets2 = union(offsets',offsets(:,1)','rows')';
end

% ---------------------------------------

function data = TwoHessian_OffsetIndex(data)

sbV = data.twoHessian.sbV;
Dimension = size(sbV,1);
sbN = size(sbV,2);
SymDimension = size(sbV,3);

refOffsets2 = data.offsets2;
nRefOffsets2 = size(refOffsets2,2);

sbV=reshape(sbV,[Dimension,sbN*SymDimension]);
[found,sbVIndex] = ismember(sbV',refOffsets2','rows');
assert(all(found));
sbVIndex = reshape(sbVIndex,[1,sbN,SymDimension]);

% Enumerate all subfacets
facetsIndex = zeros([1,sbN,SymDimension,2^SymDimension-1]);
for i=1:(2^SymDimension-1)
    b = double(dec2bin(i,SymDimension))-double('0');
    facetsIndex(1,:,:,i)=sbVIndex(1,:,b);
end
facetsIndex = reshape(facetsIndex,


data.twoHessian.sbVIndex = sbVIndex;

rmfield(data.twoHessian, 'reqOffsets2');
end

% ---------------------------------------

function data = TwoHessian_Value(data)

facetsIndex = data.twoHessian.facetsIndex;
rho = data.twoHessian.rho;
diff2=data.diff2;

Dimension=size(sbVIndex,1);
SymDimension=Dimension*(Dimension+1)/2;
sbN=size(sbVIndex,2);
nPoints=size(rho,2);

nFacets = size(facetsIndex,2);
for i=1:nFacets
    
end

end

function [value,xi,s] = TwoHessian_ElemValue(v,delta,rho)
% Minimizes <delta,xi> - rho sqrt(q0(xi)) subject to <t,xi> = 1,
% where xi represents a symmetric matrix sum_k xi_k v_k x v_k,
% tr is the trace, and q0 is dual to the two hessian.

% Output : value of the minimum, of the minimizers xi, and of sqrt(q0(xi)).
% if the minimum is not attained for xi>=0, an infinite value is returned.

    Dimension = size(v,1);
    OptDimension = size(v,2);
    nPoints = size(delta,2);
    assert(all(size(v)==[Dimension,OptDimension]));
    assert(all(size(delta)==[1,nPoints,OptDimension]));
    assert(all(size(rho)==[1,nPoints]));
    
    assert(all(rho>=0));
    
    t = zeros([OptDimension]); %trace
    q = zeros([OptDimension,OptDimension]); %two hessian
    for i=1:OptDimension
        vi=v(:,i);
        tr(i) = vi'*vi;
        for j=1:OptDimension
            vj=v(:,j);
            q(i,j) = (vi'*vi)*(vj'*vj)-(vi'*vj)^2;
        end
    end
    qR = repmat(q,[1,nPoints]);
    tR = repmat(t,[1,nPoint]);
    
    
    %Solve quadratic system
    a = t'*q*t;
    b = -VecSymVec_ScalarProduct(delta,qR,tR);
    c = VecSymVec_ScalarProduct(delta,qR,delta) - rho;
    
    disc = b.*b-a*c;
    pos = disc>=0;
    
    value = ones([1,nPoints])*Inf;
    value(pos) = (-b(pos) + sqrt(disc(pos)))./a; %TODO : which root to choose here ? Should I check both ?
    
    % Get the gradient
    dmvt = delta - delta - ScalVec_Product(values,tR);
    xi = SymVec_Product(Q, dmvt);
    xit = VecVec_ScalarProduct(xi,tR);
    xi = ScalVec_Product(1./xit, xi);

    s = ScalVec_Product(1./rho,VecVec_ScalarProduct(xi,dmvt));
    pos = pos | xi*t <=0;
    value(~pos)=Inf;
end
