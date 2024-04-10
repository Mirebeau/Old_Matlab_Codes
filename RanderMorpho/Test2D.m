addpath('/Users/mirebeau/bin/2017/ITKFM_Matlab/MFiles')

% TODO
% - Remove the Hi there (2).
% - Print the warning sign correctly. (Undisplayed error message: Seeds has inappropriate shape.)

% Create a synthetic vector field, pushing towards the unit circle.
n=100;
lin=linspace(-1.3,1.3,n);
[x,y]=meshgrid(lin,lin);
r=sqrt(x.^2+y.^2);
%imagesc(r); pause;

scale = -0.7*tanh(2*(1-r));
vx = scale.*x./r;
vy = scale.*y./r;
%imagesc(scale); pause;
imagesc(vx.*vx+vy.*vy); pause;

% Create the Rander metric.
Speed = ones(5,n,n);
Speed(1,:,:)=1;
Speed(2,:,:)=0;
Speed(3,:,:)=1;
Speed(4,:,:)=vx;
Speed(5,:,:)=vy;

clear('input');
clear('output');
input.NormType = 'Finsler2DNorm';
input.Origin=[x(1,1);y(1,1)];
input.Spacing=[x(2,2)-x(1,1); y(2,2)-y(1,1)];
input.TransposeFirstTwoImageCoordinates=1;
input.Metric=Speed;
input.Seeds = [0;0;0];
%input.Tips = [1.2,-0.5;1.2,-0.2];


output = AnisotropicFastMarching_EarlyAbort(input);

Distance = output.Distance;
%NotReached = Distance == max(Distance(:));
%Distance(NotReached) = 0;
imagesc(Distance);

% Now : extract level set, initialize and flow back.
% That will be many seeds.
