% Copyright Jean-Marie Mirebeau, 2015. jm(dot)mirebeau(at)gmail(dot)com
normBall = 1;
normValues = 1;
spacing=0.01;

[x,y]=meshgrid(-1:spacing:1, -1:spacing:1);
u = ones([1,size(x)]); u(1,:,:)=x; u(2,:,:)=y;

% ------- Direct definition of the FinslerNorm sqrt(u.m.u)-w.u -------
m = [2;0;2]; w = [0.5;0]; norm = [m;w];
if normBall; imagesc( EvaluateFinslerNorm(u,norm) <1); axis equal; pause; end;
if normValues; imagesc( EvaluateFinslerNorm(u,norm) ); axis equal; pause; end;

% --- Definition via forward speed v, reverse ratio r, side ratio s ---
v = [0.5;0.5];
r= 0.2; s=0.8;  %reverse and side speed ratios
norm = FinslerConstructor(v,r,s);
if normBall; imagesc( EvaluateFinslerNorm(u,norm) <1); axis equal; pause; end;
if normValues; imagesc( EvaluateFinslerNorm(u,norm) ); axis equal; pause; end;

% --- Likewise, but with a small side ratio s -------
v = [0.5;0.5]; 
r= 0.2; s=0.2;  %reverse and side speed ratios
norm = FinslerConstructor(v,r,s);
if normBall; imagesc( EvaluateFinslerNorm(u,norm) <1); axis equal; pause; end;
if normValues; imagesc( EvaluateFinslerNorm(u,norm) ); axis equal; pause; end;

% ---- Definition by prescribing an ellipse for the unit ball l ----
m=[4;0;1]; w=[0.2;0.3]; 
norm = TranslatedEllipse(m,w);
if normBall; imagesc( EvaluateFinslerNorm(u,norm) <1); axis equal; pause; end;
if normValues; imagesc( EvaluateFinslerNorm(u,norm) ); axis equal; pause; end;


% ---- Checking with Numerical Fast Marching -------
% Note that geodesics are straight lines, since the metric is constant.
clear('input'); clear('output');
input.NormType = 'Finsler2DNorm';
input.Origin=[x(1,1);y(1,1)];
input.Spacing=[x(2,2)-x(1,1); y(2,2)-y(1,1)];
input.TransposeFirstTwoImageCoordinates=1;

Speed = ones([5,size(x)]);
for i=1:5; Speed(i,:,:)=norm(i); end;
%  --------- !! Important !! ------------
% Fast Marching computes the distance from Tip to Seed
% In order to get the distance from Seed to Tip, one needs to reverse
% the asymmetric part w = -w; Equvalently Speed(4:5,:,:) = -Speed(4:5,:,:);
Speed(4:5,:,:) = -Speed(4:5,:,:);
input.Metric = Speed;

input.Seeds = [0;0;0];
input.Tips = [0.8,-0.5;0.8,-0.2];
output = AnisotropicFastMarching_EarlyAbort(input);
imagesc(output.Distance);
axis equal;
for i=1:size(output.Geodesics,1)
    rescaledGeodesic = RescaledCoords(output.Geodesics{i},input.Origin,input.Spacing);
    line(rescaledGeodesic(1,:),rescaledGeodesic(2,:));
end; 
