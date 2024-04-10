% Copyright Jean-Marie Mirebeau, 2015. jm(dot)mirebeau(at)gmail(dot)com

% This file illustrates the discretization of anisotropic diffusion 
% using the monotony preserving scheme published in :

% J. Fehrenbach, J.-M. Mirebeau, Sparse non-negative stencils for anisotropic diffusion,
% J. Math. Imag. Vis., vol. 49(1) (2014), pp. 123-147

% Here we diffuse a tensor field using itself.

% Important note : Matlab uses an image-as-matrix convention,
% which is not much compatible with anisotropic PDEs.
% Hence this code contains a lot of transpositions.

n=99;%n must be odd, or zero divide in this example.
[x,y]=meshgrid(-1:2/n:1,-1:2/n:1);
x=x'; y=y'; % !! transpose !!
s = size(x);

% -------- Generate a field of diffusion tensors --------
eVal1 = ones(s); eVal2 = 20*ones(s); %eigenvalues

eVec1 = zeros([s,2]); eVec2=eVec1; %eigenvectors
r=sqrt(x.*x+y.*y); 
eVec1(:,:,1)= x./r; eVec1(:,:,2)=y./r;
eVec2(:,:,1)=-y./r; eVec2(:,:,2)=x./r;

%diffusion tensors with the above eigenvectors and eigenvalues.
%Here : diffusion in a circular fashion, around the image center.
% Format : xx, xy, yy
tensors = ones([s,3]);
tensors(:,:,1)=eVal1.*eVec1(:,:,1).*eVec1(:,:,1)+eVal2.*eVec2(:,:,1).*eVec2(:,:,1);
tensors(:,:,2)=eVal1.*eVec1(:,:,1).*eVec1(:,:,2)+eVal2.*eVec2(:,:,1).*eVec2(:,:,2);
tensors(:,:,3)=eVal1.*eVec1(:,:,2).*eVec1(:,:,2)+eVal2.*eVec2(:,:,2).*eVec2(:,:,2);

%tensors(:,:,1)=1; tensors(:,:,2)=0; tensors(:,:,3)=1; %Uniform diffusion
%tensors(:,:,1)=1*(x<=0)+2*(x>0); tensors(:,:,2)=0; tensors(:,:,3)=0.01*(y<=0)+1*(y>0); %diagonal tensors

% We set identity tensors in the masked zone, 
% and anisotropic tensors in the remaining zone.

mask = (2*abs(x))<abs(y) | (2*abs(y))<abs(x) | r<0.6 | r>0.8;

tensors(:,:,1) = mask + tensors(:,:,1).*(1-mask);
tensors(:,:,2) =        tensors(:,:,2).*(1-mask);
tensors(:,:,3) = mask + tensors(:,:,3).*(1-mask);
%Showing the determinants
imagesc(tensors(:,:,1).*tensors(:,:,3)-tensors(:,:,2).^2); pause; 

% Now is the main loop.
nMainLoops=40;
nInnerLoops=10;
for mainCounter=1:nMainLoops
    % Some isotropic diffusion
    tensorsConvolved = tensors;
    for i=1:nInnerLoops
        tensorsConvolved = imfilter(tensorsConvolved,fspecial('gaussian',3,1.),'symmetric');
    end
    imagesc(tensorsConvolved(:,:,1).*tensorsConvolved(:,:,3)-tensorsConvolved(:,:,2).^2);pause;

    %Some anisotropic diffusion
    A=DiffusionSparseMatrix_2D(tensorsConvolved);
    maxTimeStep = 1./full(max(A(:)));
    dt=0.5*maxTimeStep;
    for j=1:3
        image=reshape(tensors(:,:,j)',[prod(s),1]); % !! transpose !!
        for i=1:nInnerLoops
            image=image-dt*A*image;
        end
        tensors(:,:,j)=reshape(image,s)';
    end
    imagesc(tensors(:,:,1).*tensors(:,:,3)-tensors(:,:,2).^2);
end
