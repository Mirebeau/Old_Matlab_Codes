% Copyright Jean-Marie Mirebeau, 2015. jm(dot)mirebeau(at)gmail(dot)com

% Input : 
% - two paths, approximating the right and left boundaries of an object
%    intial and final path points must coincide.
% - rho, density integrated in domain
% - sigma, weight on the curves

% options include
% - alpha in [0,1], how much of the other path is used for weights for rot solve

function [path2,path1,minVal] = ImprovePaths(path1,path2,rho,sigma,options)
    assert(all(size(rho)==size(sigma)));
    options.Spacing = [options.GridScale;options.GridScale];
    options.Size = size(rho);

    % Compute distances from paths, for weighting curl PDE
    path = [path1,fliplr(path2)];
    pathlen = size(path,2);
    dpath1 = DistanceFromPath(path(:,1:ceil(pathlen*2./3.)),      options);
    dpath2 = DistanceFromPath(path(:,floor(pathlen/3.):pathlen),  options);
    
    rho1 = (1e-6+(dpath1<options.PathWidth))./(2*options.GridScale + dpath1);
    rho2 = (1e-6+(dpath2<options.PathWidth))./(2*options.GridScale + dpath2);
    %rho1 = 1./(2*options.GridScale + dpath1); % Sad that it does not work
    %rho2 = 1./(2*options.GridScale + dpath2);
%     rho1 = 1e-6 + (dpath1<options.PathWidth);
%     rho2 = 1e-6 + (dpath2<options.PathWidth);
    
    % Solve curl PDE
    [u1,u2,p] = WeightedCurlSolve(rho,rho1,rho2,options.GridScale);
    u1SqNorm = shiftdim( u1(1,:,:).*u1(1,:,:)+u1(2,:,:).*u1(2,:,:) , 1);
    u2SqNorm = shiftdim( u2(1,:,:).*u2(1,:,:)+u2(2,:,:).*u2(2,:,:) , 1);    
    
    if GetOptions(options,'ShowDensities',false)
        imagesc(rho); title('rho'); colorbar; pause;
        imagesc(rho1); title('rho1'); colorbar; pause;
        imagesc(rho2); title('rho2'); colorbar; pause;
    end
    
    if GetOptions(options,'ShowVectorFields',false)
        quiver(options.x,options.y, shiftdim(u1(1,:,:),1) ,shiftdim(u1(2,:,:),1) ); title('u1'); pause;
        quiver(options.x,options.y, shiftdim(u2(1,:,:),1) ,shiftdim(u2(2,:,:),1) ); title('u2'); pause;
    end;
    
    if GetOptions(options,'ShowVectorFieldsNorms',false)
        imagesc(sqrt(u1SqNorm)); title('u1Norm'); colorbar; pause;
        imagesc(u1SqNorm<1); pause;
        imagesc(sqrt(u2SqNorm)); title('u2Norm'); colorbar; pause;
        imagesc(u2SqNorm<1); pause;
        imagesc(p); title('p'); colorbar;  pause;
    end
    
    %Setup Finsler metrics 
    
    Metric1 = zeros([5,options.Size]);
    Metric1(1,:,:) = max(sigma.^2, 2*u1SqNorm);
    Metric1(2,:,:)=0.;
    Metric1(3,:,:) = Metric1(1,:,:);
    Metric1(4,:,:) = u1(1,:,:); 
    Metric1(5,:,:) = u1(2,:,:);

    Metric2 = zeros([5,options.Size]);
    Metric2(1,:,:) = max(sigma.^2, 2*u2SqNorm);
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
        imagesc(d1); title('d1'); line(rescaledPath1(1,:),rescaledPath1(2,:));  pause;
        imagesc(d2); title('d2'); line(rescaledPath2(1,:),rescaledPath2(2,:));  pause;
    end
    
    if GetOptions(options,'ShowImprovedPaths',true)
        imagesc(X); title('X'); line(rescaledPath1(1,:),rescaledPath1(2,:)); line(rescaledPath2(1,:),rescaledPath2(2,:)); colorbar; pause;
    end
    
end


function d = DistanceFromPath(path,options)
    clear('input');
    input.Origin = options.Origin;
    input.Spacing = options.Spacing;
    input.NormType = 'Riemannian2DNorm';
    input.Seeds = [path; zeros(1,size(path,2))];
    input.Metric = ones([3,options.Size]);
    input.Metric(2,:,:)=0.;
    
    output = AnisotropicFastMarching_EarlyAbort(input);
    d=output.Distance;
%    imagesc(d); pause;
end