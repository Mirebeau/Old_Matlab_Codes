% Copyright Jean-Marie Mirebeau, 2015. jm(dot)mirebeau(at)gmail(dot)com

% This file is a parallel implementation of Selling's algorithm,
% for computing obtuse superbases of 2x2 SDP matrices.
% Key ingredient of the anisotropic diffusion scheme described in

% J. Fehrenbach, J.-M. Mirebeau, Sparse non-negative stencils for anisotropic diffusion,
% J. Math. Imag. Vis., vol. 49(1) (2014), pp. 123-147


function [b,weights]=TensorDecomposition_2D(tensors)

% Check data, reshape tensors
s=size(tensors);
assert(length(s)==2)
assert(s(2)==3)

% Initialize superbases
s(2)=2;
b={ones(s),ones(s),-ones(s)};
b{1}(:,2)=0;
b{2}(:,1)=0;

% Selling's algorithm
posScal=0; %number of consecutive scalar products seen
counter=0;
maxIter=50; %Could be an option.
while counter<maxIter && posScal<3
    scal = ScalarProduct_2D(b{1},b{2},tensors);
    scal = (scal>0);
    if max(scal(:))  % Superbase update: b1, -b2, b2-b1
      b{3}(scal,:) = b{2}(scal,:)-b{1}(scal,:);
      b{2}(scal,:) = -b{2}(scal,:);
      posScal=0;
    else
      posScal = posScal+1;
    end
    counter=counter+1;
    b={b{2},b{3},b{1}};
end
    
if posScal<3
    'Warning : Selling-s algorithm unterminated'
end

weights={
   -ScalarProduct_2D(b{2},b{3},tensors)/2,
   -ScalarProduct_2D(b{3},b{1},tensors)/2,
   -ScalarProduct_2D(b{1},b{2},tensors)/2
};

for i=1:3 %Compute orthogonal vectors
    c=b{i};
    b{i}(:,1)=-c(:,2);
    b{i}(:,2)= c(:,1);
end

end

function scal = ScalarProduct_2D(u,v,D)
      scal = D(:,1).*u(:,1).*v(:,1) +  D(:,2).*(u(:,1).*v(:,2) + u(:,2).*v(:,1) ) + D(:,3).*u(:,2).*v(:,2);
end