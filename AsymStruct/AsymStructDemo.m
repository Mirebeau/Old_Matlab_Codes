% Copyright Jean-Marie Mirebeau, 2015. jm(dot)mirebeau(at)gmail(dot)com

addpath('/Users/mirebeau/bin/2015/Matlab/MFiles') % -> AnisotropicFastMarching_EarlyAbort.mexmaci64

addpath('../Eig3Folder')
addpath('../GetOptions')
addpath('../ParallelMatrix')
n=100;
interval = linspace(-1.2,1.2,n);
[x,y]=meshgrid(interval,interval);

%Image is noisy characteristic function of the unit ball. 
image = x.^2+y.^2<1;
imageSize = size(image);
noiseLevel=0.1;
image = image+noiseLevel*rand(imageSize);
clf;
imagesc(image);
saveas(gcf,'Output/Image','png');

% Build Finsler norm
sdpExponent = 1; % Choose <1 to reduce anisotropy.
clear('options');
options.rotate=true;
options.eps_absolute = 0.02;
finsler = AsymStruct(image,options);
finsler = reshape(finsler,[5,prod(imageSize)]);
sdp = FinslerFromToSDP(finsler);
sdp = SymmetricMatrix3Power(sdp, -sdpExponent); %minus is for duality %Safe option needed for identities...
finsler = FinslerFromToSDP(sdp);
%finsler(1,:)=1;finsler(2,:)=0;finsler(3,:)=1;finsler(4,:)=0;finsler(5,:)=0; %Set Identity for testing
finsler=reshape(finsler,[5,imageSize]);
dlmwrite('Output/FinslerNorms.txt',Transpose_3D(finsler));

% setup fast marching
input.Origin=[x(1,1);y(1,1)];
input.Spacing=[x(2,2)-x(1,1); y(2,2)-y(1,1)];
input.TransposeFirstTwoImageCoordinates=1; %!!
input.NormType = 'Finsler2DNorm';
input.Metric = finsler;
input.Seeds = [-1;0;0];
input.Tips = [0;-1];
%input.Tips = [0,0;1,-1];
input.StopWhenTipsAreReached = 1;

output = AnisotropicFastMarching_EarlyAbort(input);

% To get a better image, set to zero the points not reached by the Fast
% Marching algorithm (by default DBL_MAX)
Distance = output.Distance;
NotReached = Distance == max(Distance(:));
Distance(NotReached) = 0;

%Display
pause;
clf;
imagesc(Distance);
saveas(gcf,'Output/Distance','png');

for i=1:size(output.Geodesics,1)
    rescaledGeodesic = RescaledCoords(output.Geodesics{i},input.Origin,input.Spacing);
    line(rescaledGeodesic(1,:),rescaledGeodesic(2,:));
end;
fprintf('\n\n');
