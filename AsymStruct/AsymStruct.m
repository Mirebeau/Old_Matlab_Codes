% Copyright Jean-Marie Mirebeau, 2015. jm(dot)mirebeau(at)gmail(dot)com

% Asymmetric variant of the image structure tensor, 
% based on Finsler asymmetric norms.
% Options : 
% - noise_filter, feature_filter -> gaussian regularization filters, as in structure tensor
% - eps_relative, eps_absolute -> avoid tensor degeneracy by adding ...
%   (eps_relative+eps_absolute * Tr(M))*Id
% - rotate : rotate gradients by pi/2, thus also norms, for applications to image segmentation.

% - TransposeGradient : related to Matlab's image as matrix convention,
% and to similar option in AnisotropicFastMarching_EarlyAbort 
% Compatibility requires TransposeGradient = ~ TransposeFirstTwoImageCoordinates

function finsler = AsymStruct(image,options)
    %convolve at the noisescale
    image = imfilter(image,...
      GetOptions(options,'noise_filter',fspecial('gaussian',3,0.7)),...
      'symmetric');
  
    if(GetOptions(options,'TransposeGradient',0))
      [gx,gy]=gradient(image'); gx=gx'; gy=gy';
    else
      [gx,gy]=gradient(image); 
    end
    
    if(GetOptions(options,'rotate',0))
       temp = gx;
       gx=-gy;
       gy=temp;
    end
    
    %Regularization via 'articificial noise'
    eps = (gx.*gx+gy.*gy)*GetOptions(options,'eps_relative',0.)...
        +GetOptions(options,'eps_absolute',0.);

    %Assembling the metric
    finsler = ones([3,size(image)]);
    finsler(1,:,:)=gx.*gx + eps;
    finsler(2,:,:)=gx.*gy;
    finsler(3,:,:)=gy.*gy + eps;
    finsler(4,:,:)=gx;
    finsler(5,:,:)=gy;
    
    %Regularization at feature scale
    featureFilter =  GetOptions(options,'feature_filter',fspecial('gaussian',5,1.5));
    for i=1:5 finsler(i,:,:)=imfilter(finsler(i,:,:),featureFilter,'symmetric'); end;
  
end