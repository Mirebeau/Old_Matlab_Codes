% Copyright Jean-Marie Mirebeau, 2016. jm(dot)mirebeau(at)gmail(dot)com

% This function takes as input a Finsler metric (sqrt(quadratic)+linear form).
% Quadratic form is given by m11,m21,m22,m31,m32,m33, ... Linear w1,w2,w3,...
% Definiteness assumption : w w^T < m in the sense of matrices

% options.dims gives the domain dimensions, of length 1,2 or 3. Pixel size is 1.
% options.eps gives the resolvant scaling
% Output is a sparse matrix specification A, discretizing a second operator, 
% which inverse kernel A^-1(x,y) approximates exp(-d_F(x,y)/eps)


function [s,w,a] = EikoFinslerOperator(m,w,eps)

Dimension = size(w,1);
nPoints = size(w,2);
SymDimension = Dimension*(Dimension+1)/2;

assert(all(size(m) == [SymDimension,nPoints]));
assert(all(size(w) == [Dimension,nPoints]));

% Replace metric (meant to measure vectors) 
% with the dual finsler metric (meant to measure gradients),
% as in the eikonal equation.
% Dual metric still obeys positivity assumption.

[m,w] = SymVec_DualFinsler(m,w);

% Construct conductivity tensors s=m-w w^T, which are positive by construction

s = m - Vec_SelfOuterProduct(w);

% Recale the different constributions.

a = ones([1,nPoints]); 

w = -2*eps*w;
s = eps^2*s;

end