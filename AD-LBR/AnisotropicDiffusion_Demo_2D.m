% Copyright Jean-Marie Mirebeau, 2015. jm(dot)mirebeau(at)gmail(dot)com

% This file illustrates the discretization of anisotropic diffusion 
% using the monotony preserving scheme published in :

% J. Fehrenbach, J.-M. Mirebeau, Sparse non-negative stencils for anisotropic diffusion,
% J. Math. Imag. Vis., vol. 49(1) (2014), pp. 123-147

% Diffusion tensors are here built by hand, instead of following
% J. Weickert's structure tensor based constructions.

% Important note : Matlab uses an image-as-matrix convention,
% which is not much compatible with anisotropic PDEs.
% Hence this code contains a lot of transpositions.

n=99;%n must be odd, or zero divide.
[x,y]=meshgrid(-1:2/n:1,-1:2/n:1);
x=x'; y=y'; % !! transpose !!
s = size(x);

% -------- Generate a field of diffusion tensors --------
eVal1 = 0.05*ones(s); eVal2 = 1*ones(s); %eigenvalues

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

% ---------- Generate the operator matrix ----------
A=DiffusionSparseMatrix_2D(tensors);
maxTimeStep = 1./full(max(A(:)));
dt=0.5*maxTimeStep;

disp('----- Demo : anisotropic diffusion of some noise, circularly. ------')
image = random('Uniform',0,1,s(2),s(1)); % diffused image, matlab convention

%Matrix vector product requires reshape
image=reshape(image',[prod(s),1]); % !! transpose !!
for i=1:10
    image=image-dt*A*image;
end
image=reshape(image,s)'; % !! transpose !!

imshow(image)
pause()

disp('----- Demo : distance map from heat, w.r.t to an anisotropic metric, using Varadhan formula. ------');
seed=[10,40]; % seed for distance computations. (xCoord, yCoord), in pixels.
image=0.*image;
image(seed(2),seed(1))=1;

%Matrix vector product requires reshape
image=reshape(image',[prod(s),1]); % !! transpose !!
eps = 0.01; 
image = (speye(prod(s)) + eps*A)\image; 
image=reshape(image,s)'; % !! transpose !!

image=-log(image);
image=image-min(image(:));
image=image/max(image(:));

image(seed(2),seed(1))=1; 
% seed is a white dot.
% black : close to seed, for the riemannian distance.
% white : far from seed, for the riemannian distance.
imshow(image)