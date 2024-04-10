% Copyright Jean-Marie Mirebeau, 2016. jm(dot)mirebeau(at)gmail(dot)com

% This function takes as input 
% - s : positive definite tensor field
% - w : vector field
% - a : scalar field
% - options, with fields specifying : 
%       * dims : the grid dimensions.
%       * gridscale : the pixel scale.
%       * nonDivergence : whether to use the divergence form for the second order part. (default false)
%       * dirichlet : wether to use null Dirichlet boundary conditions, Neumann otherwise. (default false)

% It assembles the finite difference matrix of 
% a u + <w,grad(u)> - Tr(s Hessian(u))
% or if the divergence form is required
% a u + <w,grad(u)> - div(s grad(u))

% The second order operator discretization is based on the AD-LBR scheme
% J. Fehrenbach, J.-M. Mirebeau, Sparse non-negative stencils for anisotropic diffusion,
% J. Math. Imag. Vis., vol. 49(1) (2014), pp. 123-147

% The first order operator discretization is upwind, and uses the same stencils.

% (Caution : Not too much attention is paid to the boundary conditions.
% In particular, not very sure that Neumann with nonDivergence has much meaning, or for the first order part.
% For the geodesics in heat application, the function decays rapidly, so this should not be an issue.) 

function [Rows,Cols,Vals] = DifferentialOperatorMatrix(s,w,a,options)

dims = options.dims;
Dimension = size(dims,2);
SymDimension = Dimension*(Dimension+1)/2;
nPoints = prod(dims);
h = options.gridScale;
dirichlet = isfield(options,'dirichlet') && options.dirichlet;
nonDivergence = isfield(options,'nonDivergence') && options.nonDivergence;

% ------- Build the linear operator ------

% Construct the stencils
[offsets,weights] = Sym_Decomposition(s);
nOffsets = size(offsets,3);
indices = zeros([1,nPoints,nOffsets,2]);
indices(:,:,:,1) = GridIndices(dims, offsets);
indices(:,:,:,2) = GridIndices(dims,-offsets);

Rows = [];
Cols = [];
Vals = [];

if numel(a)>0
    %Order zero : a*u
    assert(all(size(a)==[1,nPoints]));
    Rows = [Rows, 1:nPoints];
    Cols = [Cols, 1:nPoints];
    Vals = [Vals, a];
end

if numel(w)>0
    %Order one : 2 eps <w, grad u> = 2 eps <v,s grad(u)>
    assert(all(size(w)==[Dimension,nPoints]));
    v = SymVec_Solve(s,w);

    for iOffset = 1:nOffsets
        coef = h^(-1) * VecVec_ScalarProduct(v,offsets(:,:,iOffset)) .* weights(1,:,iOffset);
        
        neigh = zeros([1,nPoints]);
        neigh(coef>=0) = indices(1,coef>=0,iOffset,1);
        neigh(coef<0)  = indices(1,coef<0, iOffset,2);
        
        coef = abs(coef);
        pos = 1:nPoints;
        if ~dirichlet
            inbox   = (neigh ~= 0);
            pos     = pos(inbox);
            neigh   = neigh(inbox);
            coef    = coef(inbox);
        end
        Rows = [Rows,   pos,    pos];
        Cols = [Cols,   pos,    neigh];
        Vals = [Vals,   coef,   -coef]; 
    end
end

if numel(s)>0
    %Order two : -Tr(s Hessian(u)) or -eps^2 div( s grad(u)), which is the matrix of <grad(u), s grad(u)>
    assert(all(size(s)==[SymDimension,nPoints]));
    for iOffset = 1:nOffsets
        for i=1:2
            coef = h^(-2) * weights(1,:,iOffset);
            neigh = indices(1,:,iOffset,i);
            pos = 1:nPoints;
                        
            if ~dirichlet
                inbox   = (neigh~=0);
                pos     = pos(inbox);
                neigh   = neigh(inbox);
                coef    = coef(inbox);
            end
            
            if ~nonDivergence
                Rows = [Rows,   pos,    pos,    neigh,  neigh];
                Cols = [Cols,   pos,    neigh,  pos,    neigh];
                Vals = [Vals,   coef/2, -coef/2,-coef/2,coef/2];
            else
                Rows = [Rows,   pos,    pos];
                Cols = [Cols,   pos,    neigh];
                Vals = [Vals,   coef,   -coef];
            end
        end
    end
end

inbox   = Rows~=0 & Cols~=0 & Vals~=0; %Out of box index or null matrix coefficient.
Rows = Rows(inbox);
Cols = Cols(inbox);
Vals = Vals(inbox);

end