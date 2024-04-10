% Copyright Jean-Marie Mirebeau, 2015. jm(dot)mirebeau(at)gmail(dot)com

% Input : 
% - two paths, approximating the right and left boundaries of an object
%    intial and final path points must coincide.
% - rho, density integrated in domain
% - sigma, weight on the curves

% options include
% ??- alpha in [0,1], how much of the other path is used for weights for rot solve
% - ShowDistanceFromPath, ShowDensities, ShowVectorFields (costly), ShowVectorFieldsNorms, ShowDistanceMaps, ShowImprovedPaths

function [path2,path1,minVal] = ImprovePaths_loc(path1,path2,rho,sigma,options)
    assert(all(size(rho)==size(sigma)));
    options.Spacing = [options.GridScale;options.GridScale];
    options.Size = size(rho);

    % Compute distances from paths
    distWeight = abs(rho)./sigma;
    distWeight = distWeight + mean(distWeight(:))/100;
    path = [path1,fliplr(path2)];
    pathlen = size(path,2);
    n1 = floor(pathlen/3.); 
    n2 = ceil(pathlen*2./3.);
    
    neighborhoodWidth=0.5;
    dpath1 = DistanceFromPath(path(:,1:n2),       distWeight,  options, neighborhoodWidth);
    dpath2 = DistanceFromPath(path(:,n1:pathlen), distWeight,  options, neighborhoodWidth);
    dpath12 = DistanceFromPath(path(:,n1:n2),     distWeight,  options, neighborhoodWidth);
    
    if GetOptions(options,'ShowDistanceFromPath',false)
        farValue = realmax('double')/2.;
        dp=dpath1;  dp(dp==farValue)=0; imagesc(dp); title('dpath1');  pause;
        dp=dpath2;  dp(dp==farValue)=0; imagesc(dp); title('dpath2');  pause;
        dp=dpath12; dp(dp==farValue)=0; imagesc(dp); title('dpath12'); pause;
    end;

    % Setup weights for curl PDE
    
    rho1 = dpath1 < neighborhoodWidth; 
    rho2 = dpath2 < neighborhoodWidth; 
    rhop = dpath12 < neighborhoodWidth; rhop=rhop.*rho1.*rho2; %Only positivity matters
    
    rho1=rho1./(2*options.GridScale+dpath1);
    rho2=rho2./(2*options.GridScale+dpath2);
    
%    rho1 = 1./(2*options.GridScale + dpath1);
%    rho2 = 1./(2*options.GridScale + dpath2);
%    rho1 = 1e-6 + (dpath1<options.PathWidth);
%    rho2 = 1e-6 + (dpath2<options.PathWidth);

    
    if GetOptions(options,'ShowDensities',false)
        imagesc(rho); title('rho'); colorbar; pause;
        imagesc(rho1); title('rho1'); colorbar; pause;
        imagesc(rho2); title('rho2'); colorbar; pause;
        imagesc(rhop); title('rhop'); colorbar; pause;
    end;

    
    % Solve curl PDE
    [u1,u2,p] = WeightedCurlSolve_Loc(rho,rho1,rho2,options.GridScale,rhop);
    u1SqNorm = shiftdim( u1(1,:,:).*u1(1,:,:)+u1(2,:,:).*u1(2,:,:) , 1);
    u2SqNorm = shiftdim( u2(1,:,:).*u2(1,:,:)+u2(2,:,:).*u2(2,:,:) , 1);    
    
    
    if GetOptions(options,'ShowVectorFields',false)
        quiver(options.x,options.y, shiftdim(u1(1,:,:),1) ,shiftdim(u1(2,:,:),1) ); title('u1'); pause;
        quiver(options.x,options.y, shiftdim(u2(1,:,:),1) ,shiftdim(u2(2,:,:),1) ); title('u2'); pause;
    end;
    
    if GetOptions(options,'ShowVectorFieldsNorms',false)
        imagesc(sqrt(u1SqNorm)); title('u1Norm'); colorbar; pause;
        imagesc(u1SqNorm<1); title('u1Norm<1'); pause;
        imagesc(sqrt(u2SqNorm)); title('u2Norm'); colorbar; pause;
        imagesc(u2SqNorm<1); title('u2Norm<1'); pause;
        imagesc(p); title('p'); colorbar;  pause;
    end;
    
    %Setup Finsler metrics 
    % TO DO : check that 2*u1SqNorm < sigma.^2 on a neighborhood of the curve, e.g. the set dpath1<(neighborhoodWidth/2). 
    % Otherwise e.g. set neighborhoodWidth=neighborhoodWidth/2; and go back to 'Setup weights for curl PDE'.
    
    Metric1 = zeros([5,options.Size]);
    Metric1(1,:,:) = max(sigma.^2, 2*u1SqNorm); 
    domain1=1./rho1; domain1(domain1<Inf)=0; Metric1(1,:,:) = max(shiftdim(Metric1(1,:,:),1),domain1); %Restrict FM to domain1
    Metric1(2,:,:)=0.;
    Metric1(3,:,:) = Metric1(1,:,:);
    Metric1(4,:,:) = u1(1,:,:); 
    Metric1(5,:,:) = u1(2,:,:);

    Metric2 = zeros([5,options.Size]);
    Metric2(1,:,:) = max(sigma.^2, 2*u2SqNorm); 
    domain2=1./rho2; domain2(domain2<Inf)=0; Metric2(1,:,:) = max(shiftdim(Metric2(1,:,:),1),domain2);
    Metric2(2,:,:)=0.;
    Metric2(3,:,:) = Metric2(1,:,:);
    Metric2(4,:,:) = u2(1,:,:);
    Metric2(5,:,:) = u2(2,:,:);

    
    % Setup and run fast marching
    clear('input');
    input.Origin = options.Origin;
    input.Spacing=options.Spacing;
    input.NormType = 'Finsler2DNorm';
    input.Seeds = [path(:,1);0];
    
    input.Metric = Metric1;
    output = AnisotropicFastMarching_EarlyAbort(input);
    d1 = output.Distance;
    
    input.Metric = Metric2;
    output = AnisotropicFastMarching_EarlyAbort(input);
    d2 = output.Distance;

    % Find minimum value
    X = d1+d2-p;
    X(rhop==0)=Inf;
    [minVal,ind] = min(X(:));
    [I,J] = ind2sub(options.Size,ind);
    tip = RescaledCoords([J;I],input.Origin,input.Spacing,false);
    tip
    X(I,J) = max(X(:));
    
    % Re-run fast marching to get the paths
    input.StopWhenTipsAreReached = 1;
    input.Tips = tip;
    
    input.Metric = Metric1;
    output = AnisotropicFastMarching_EarlyAbort(input);
    path1 = output.Geodesics{1};

    input.Metric = Metric2;
    output = AnisotropicFastMarching_EarlyAbort(input);
    path2 = output.Geodesics{1};
    

    rescaledPath1 = RescaledCoords(path1,input.Origin,input.Spacing);
    rescaledPath2 = RescaledCoords(path2,input.Origin,input.Spacing);
    
    if GetOptions(options,'ShowDistanceMaps',false)
        dp=d1; dp(dp==farValue)=0; imagesc(dp); title('d1'); line(rescaledPath1(1,:),rescaledPath1(2,:));  pause;
        dp=d2; dp(dp==farValue)=0; imagesc(dp); title('d2'); line(rescaledPath2(1,:),rescaledPath2(2,:));  pause;
    end;
    
    if GetOptions(options,'ShowImprovedPaths',true)
        imagesc(X); title('X'); line(rescaledPath1(1,:),rescaledPath1(2,:)); line(rescaledPath2(1,:),rescaledPath2(2,:)); colorbar; pause;
    end;
    
end


function d = DistanceFromPath(path,weight,options, neighborhoodWidth)
    clear('input');
    input.Origin = options.Origin;
    input.Spacing = options.Spacing;
    input.NormType = 'Riemannian2DNorm';
    input.Seeds = [path; zeros(1,size(path,2))];
    input.Metric = zeros([3,options.Size]);
    input.Metric(1,:,:)=weight.^2;
    input.Metric(3,:,:)=input.Metric(1,:,:);
    input.StopAtDistance = neighborhoodWidth;
    
    output = AnisotropicFastMarching_EarlyAbort(input);
    d=output.Distance;
%    imagesc(d); pause;
end