% Copyright Jean-Marie Mirebeau, 2015. jm(dot)mirebeau(at)gmail(dot)com
% This file implements in 2D anisotropic divergence form diffusion, 
% with prescribed diffusion tensors, using the scheme described in : 

% J. Fehrenbach, J.-M. Mirebeau, Sparse non-negative stencils for anisotropic diffusion,
% J. Math. Imag. Vis., vol. 49(1) (2014), pp. 123-147

function A=DiffusionSparseMatrix_2D(tensors)

  %first compute the obtuse superbases
  s = size(tensors);
  s = s([1,2]); %image size
  [superbases,weights] = TensorDecomposition_2D(reshape(tensors, [prod(s),3]));
  
  %Make the sparse matrix
  indices_range = 1:prod(s);
  column=[];
  row=[];
  coefficient=[];
  
  [x,y] = meshgrid( 1:s(1), 1:s(2) ); 
  x=x'; y=y'; %careful with axes ...
  z=ones([s,2]);
  z(:,:,1)=x;
  z(:,:,2)=y;
  z=reshape(z,[prod(s),2]);
  
  for eps = (-1):2:1
    for i=1:3
      neighbor = z+eps*superbases{i};
            
      inside = neighbor(:,1)>=1 & neighbor(:,1)<=s(1) & neighbor(:,2)>=1 & neighbor(:,2)<=s(2);
      scal=weights{i}; 
      scal = scal'; % line to column, not an axes issue
      
      neighbor = neighbor(inside,:);
      points = indices_range(inside);
      scal = scal(inside);
      
      neighbor = neighbor(:,1)+s(1)*(neighbor(:,2)-1); %conversion to 1D index
      neighbor = neighbor'; % line to column, not an axes issue
      
      column =      [column, points];
      row =         [row, points];
      coefficient = [coefficient, scal];
        
      column      = [column,      points,    neighbor,  neighbor];
      row         = [row,         neighbor,  points,    neighbor];
      coefficient = [coefficient, -scal,     -scal,     scal];
    end 
  end
  
    A=sparse(column, row, coefficient, prod(s), prod(s) );
end