% Copyright Jean-Marie Mirebeau, 2015. jm(dot)mirebeau(at)gmail(dot)com

% This file implements 3D anisotropic divergence form diffusion, 
% with prescribed diffusion tensors, using the scheme described in : 

% J. Fehrenbach, J.-M. Mirebeau, Sparse non-negative stencils for anisotropic diffusion,
% J. Math. Imag. Vis., vol. 49(1) (2014), pp. 123-147

% Remark on performance: On 'large' cases, computation time 
% is dominated by the sparse matrix assembly "spmat(col,row,coef,n,n)". 
% In case of need, consider the following optimized C++ implementation designed for
% the Insight Toolkit (ITK) 

% J. Fehrenbach, J.-M. Mirebeau, L. Risser, S. Tobji,
% Anisotropic Diffusion in ITK, Insight Journal, 2015
% http://www.insight-journal.org/browse/publication/953

function A=DiffusionSparseMatrix_3D(tensors, periodic3)
  if nargin<2
      periodic3=false; 
  end
  %first compute the obtuse superbases
  s = size(tensors);
  s = s(1:3); %image size

  [superbases,weights] = TensorDecomposition_3D(reshape(tensors, [prod(s),6]));

  %Make the sparse matrix
  indices_range = 1:prod(s);
  column=[];
  row=[];
  coefficient=[];
  
  [x,y,z] = meshgrid( 1:s(1), 1:s(2), 1:s(3) ); 
  x=Transpose_3D(x); y=Transpose_3D(y); z=Transpose_3D(z); %careful with axes ...
  w=zeros([s,3]);
  w(:,:,:,1)=x;
  w(:,:,:,2)=y;
  w(:,:,:,3)=z;
  w=reshape(w,[prod(s),3]);

  for eps = (-1):2:1
    for i=1:6
      neighbor = w+eps*superbases{i};
            
      if(~periodic3)
        inside = ...
            neighbor(:,1)>=1 & neighbor(:,1)<=s(1) & ... 
            neighbor(:,2)>=1 & neighbor(:,2)<=s(2) & ...
            neighbor(:,3)>=1 & neighbor(:,3)<=s(3);
        neighbor = neighbor(inside,:);
      else 
        inside = ...
            neighbor(:,1)>=1 & neighbor(:,1)<=s(1) & ... 
            neighbor(:,2)>=1 & neighbor(:,2)<=s(2); %Periodicity only affects the third axis
        neighbor = neighbor(inside,:);
        neighbor(:,3)=1+mod(neighbor(:,3)-1, s(3));
      end

      scal = weights{i}; 
      scal = scal'; % line to column, not an axes issue
      
      points = indices_range(inside);
      scal = scal(inside);
      
      neighbor = neighbor(:,1)+s(1)*((neighbor(:,2)-1) + s(2)*(neighbor(:,3)-1)); %conversion to 1D index
      neighbor = neighbor'; % line to column, not an axes issue
      
      column =      [column, points];
      row =         [row, points];
      coefficient = [coefficient, scal];
        
      column      = [column,      points,    neighbor,  neighbor];
      row         = [row,         neighbor,  points,    neighbor];
      coefficient = [coefficient, -scal,     -scal,     scal];
    end 
  end
  
  % Surprisingly, the bottleneck, 70% of cpu time in large cases, is actually building the 
  % sparse matrix from the triplets. No way to escape this.
  disp('Sparse matrix construction');
  tic
  A=sparse(column, row, coefficient, prod(s), prod(s) );
  toc
end