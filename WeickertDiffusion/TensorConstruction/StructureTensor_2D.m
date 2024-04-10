% Copyright Jean-Marie Mirebeau, 2015. jm(dot)mirebeau(at)gmail(dot)com

function tensor = StructureTensor_2D(image,options)

  if ismatrix(image)
      tensor = StructureTensor_SingleChannel_2D(image,options);
  elseif ndims(image)==3
    % For multichannel images, use the sum of layer's structure tensors
    tensor = zeros([size(image,1),size(image,2),3]);
    for layer = 1:size(image,3)
        tensor = tensor+StructureTensor_SingleChannel_2D(image(:,:,layer),options);
    end
  else 
      error('Structure tensor: invalid number of image dimensions');
  end;

  if GetOptions(options,'rescale_for_unit_maximum_trace',1)
      maxTrace = max(tensor(:,1)+tensor(:,3));
      tensor = tensor/maxTrace;
  end
  
  
end

function tensor = StructureTensor_SingleChannel_2D(image,options)
  image = imfilter(image,...
      GetOptions(options,'noise_filter',fspecial('gaussian',3,0.7)),...
      'symmetric');
  
  [gx,gy]=gradient(image');
  gx=gx';
  gy=gy';
  
  tensor = ones([size(image),3]);
  tensor(:,:,1)=gx.*gx;
  tensor(:,:,2)=gx.*gy;
  tensor(:,:,3)=gy.*gy;
  tensor = imfilter(tensor,...
      GetOptions(options,'feature_filter',fspecial('gaussian',5,1.5)),...
      'symmetric'); 
end