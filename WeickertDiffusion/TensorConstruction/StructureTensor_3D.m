% Copyright Jean-Marie Mirebeau, 2015. jm(dot)mirebeau(at)gmail(dot)com

function tensor = StructureTensor_3D(image,options)

  if ndims(image)==3
      tensor = StructureTensor_SingleChannel_3D(image,options);
  elseif ndims(image)==4
    % For multichannel images, use the sum of layer's structure tensors
    tensor = zeros([size(image,1),size(image,2),size(image,3),6]);
    for layer = 1:size(image,4)
        tensor = tensor+StructureTensor_SingleChannel_3D(image(:,:,layer),options);
    end
  else 
      error('Structure tensor: invalid number of image dimensions');
  end;

  if GetOptions(options,'rescale_for_unit_maximum_trace',1)
      maxTrace = max(tensor(:,1)+tensor(:,3)+tensor(:,6));
      tensor = tensor/maxTrace;
  end
  
end

function tensor = StructureTensor_SingleChannel_3D(image,options)

  image=GaussianFilter_3D(image,...
      GetOptions(options,'noise_filter',fspecial('gaussian',[3,1],0.7)));
  
  [gx,gy,gz]=gradient(Transpose_3D(image));
  gx=Transpose_3D(gx);
  gy=Transpose_3D(gy);
  gz=Transpose_3D(gz);
  
  tensor = zeros([size(image),6]);
  tensor(:,:,:,1)=gx.*gx;
  tensor(:,:,:,2)=gx.*gy;
  tensor(:,:,:,3)=gy.*gy;
  tensor(:,:,:,4)=gx.*gz;
  tensor(:,:,:,5)=gy.*gz;
  tensor(:,:,:,6)=gz.*gz;
  tensor = GaussianFilter_3D(tensor,...
      GetOptions(options,'feature_filter',fspecial('gaussian',[5,1],1.5))); 
end