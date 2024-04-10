% Copyright Jean-Marie Mirebeau, 2015. jm(dot)mirebeau(at)gmail(dot)com

% This file is a parallel implementation of Selling's algorithm,
% for computing obtuse superbases of 3x3 SDP matrices.
% Key ingredient of the anisotropic diffusion scheme described in

% J. Fehrenbach, J.-M. Mirebeau, Sparse non-negative stencils for anisotropic diffusion,
% J. Math. Imag. Vis., vol. 49(1) (2014), pp. 123-147

function [offsets,weights]=TensorDecomposition_3D(tensors)

% Check data, reshape tensors
s=size(tensors);
assert(length(s)==2)
assert(s(2)==6)

% Initialize superbases
s(2)=3;
b={zeros(s),zeros(s),zeros(s),-ones(s)};
b{1}(:,1)=1;
b{2}(:,2)=1;
b{3}(:,3)=1;

% Selling's algorithm
counter=0;
maxIter=200; %increase for highly anisotropic tensor fields
done=0;
while done==0
    scalPos=0;
    for i=1:3
        for j=(i+1):4
            counter=counter+1;
            if counter>=maxIter,
                disp('Warning : Selling-s algorithm unterminated')
                done=1;
                break;
            end       
            
            scal = ScalarProduct(b{i},b{j},tensors);
            scal = (scal>0);
            if max(scal(:)) %Superbase update -b1, b2, b3+b1, b4+b1
                for k=1:4
                    if k~=i && k~=j
                        b{k}(scal,:)=b{k}(scal,:)+b{i}(scal,:);
                    end
                end
                b{i}(scal,:)=-b{i}(scal,:);
                scalPos=1;
            end
            if i==3 && j==4 && ~scalPos, done=1; end
        end
        if done, break, end
    end
end

weights={ ...
   -ScalarProduct(b{1},b{2},tensors)/2, ...
   -ScalarProduct(b{1},b{3},tensors)/2, ...
   -ScalarProduct(b{1},b{4},tensors)/2, ...
   -ScalarProduct(b{2},b{3},tensors)/2, ...
   -ScalarProduct(b{2},b{4},tensors)/2, ...
   -ScalarProduct(b{3},b{4},tensors)/2 ...
};

offsets = { ...
    CrossProduct(b{3},b{4}), ...
    CrossProduct(b{2},b{4}), ...
    CrossProduct(b{2},b{3}), ...
    CrossProduct(b{1},b{4}), ...
    CrossProduct(b{1},b{3}), ...
    CrossProduct(b{1},b{2}) ...
};

for i=1:3 %Compute cross products
    c=b{i};
    b{i}(:,1)=-c(:,2);
    b{i}(:,2)= c(:,1);
end

end

function cross = CrossProduct(u,v)
    cross = zeros(size(u));
    cross(:,1)=u(:,2).*v(:,3)-u(:,3).*v(:,2);
    cross(:,2)=u(:,3).*v(:,1)-u(:,1).*v(:,3);
    cross(:,3)=u(:,1).*v(:,2)-u(:,2).*v(:,1);
end

% Format : xx, xy, yy, xz, yz, zz
function scal = ScalarProduct(u,v,D)
      scal = ... %
        D(:,1).*u(:,1).*v(:,1) + ... %xx
        D(:,2).*(u(:,1).*v(:,2) + u(:,2).*v(:,1) ) + ... %xy
        D(:,3).*u(:,2).*v(:,2) + ... %yy
        D(:,4).*(u(:,1).*v(:,3) + u(:,3).*v(:,1) ) + ... %xz
        D(:,5).*(u(:,2).*v(:,3) + u(:,3).*v(:,2) ) + ... %yz
        D(:,6).*u(:,3).*v(:,3); %zz
end