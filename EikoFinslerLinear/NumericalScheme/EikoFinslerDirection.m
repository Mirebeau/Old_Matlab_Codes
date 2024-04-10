% Copyright Jean-Marie Mirebeau, 2016. jm(dot)mirebeau(at)gmail(dot)com

% This function outputs the distorted gradient direction, 
% associated to a Finsler metric, and used to compute geodesics.
% It uses the same adaptive stencils as EikoFinslerLinear

% It takes as input some ouputs of EikoFinslerLinear 
% - s : positive definite tensor field
% - w : vector field
% -sol: solution to the linear system
% -options : 
%       * dims 
%       * gridScale

function [dir] = EikoFinslerDirection(s,w,u,options)

dims = options.dims;
Dimension = size(dims,2);
SymDimension = Dimension*(Dimension+1)/2;
nPoints = prod(dims);
h = options.gridScale;
eps = options.eps;

assert(all(size(s)==[SymDimension,nPoints]));
assert(all(size(w)==[Dimension,nPoints]));
assert(all(size(u)==[1,nPoints]));

s = s/eps^2;
w = w/(2*eps);

% * offsets,weights : (?!TODO!?) precomputed Sym_Decomposition(s). (from EikoFinslerLinear, given through options.)
[offsets,weights] = Sym_Decomposition(s);
nOffsets = size(offsets,3);
indices = zeros([1,nPoints,nOffsets,2]);
indices(:,:,:,1) = GridIndices(dims, offsets);
indices(:,:,:,2) = GridIndices(dims,-offsets);

% s*grad(u) and grad(u)*s*grad(u)
sgu = zeros([Dimension,nPoints]);
n2sgu = zeros([1,nPoints]);

for iOffset=1:nOffsets
    for i=1:2
        neigh = indices(1,:,iOffset,i);
        revNeigh = indices(1,:,iOffset,3-i);
        
        inbox    = neigh~=0;
        revInbox = revNeigh~=0;
        
        coef = h^(-1) * (1-0.5*revInbox) .* weights(1,:,iOffset);        
        sgu(:,inbox) = sgu(:,inbox)+ScalVec_Product((3-2*i)*coef(inbox).*(u(neigh(inbox))-u(inbox)),offsets(:,inbox,iOffset));
        n2sgu(inbox)=n2sgu(inbox)+coef(inbox).*(u(neigh(inbox))-u(inbox)).^2;
    end
end

% *m*grad(u) = (s+ww^T)*grad u = sgu + w*<v,sgu>.
v = SymVec_Solve(s,w);
vsgu = VecVec_ScalarProduct(v,sgu);
mgu = sgu + ScalVec_Product(vsgu,w);

n2mgu = n2sgu+vsgu.^2;
dir = ScalVec_Product(1./sqrt(n2mgu),mgu) - w; % So, plus or minus ?
end