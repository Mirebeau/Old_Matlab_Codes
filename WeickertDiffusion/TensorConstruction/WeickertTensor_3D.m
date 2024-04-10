% Copyright Jean-Marie Mirebeau, 2015. jm(dot)mirebeau(at)gmail(dot)com

function wt = WeickertTensor_3D(tensor,options)

  if(ndims(tensor)~=4 || size(tensor,4)~=6) error('WeickertTensor_3D: invalid tensor format'); end;
  
  lambda = GetOptions(options,'Weickert_lambda');
  alpha = GetOptions(options,'Weickert_alpha',0.1);
  m = GetOptions(options,'Weickert_m',2);
  choice = GetOptions(options,'Weickert_choice','cEED');
  
  s=size(tensor); s=s(1:3); prods = prod(s);
  % Using Eig3 for solving multiple eigenvalues.
  % Symmetric tensors
  c = reshape(tensor,[prods,6]);
  st = zeros([3,3,prods]);
  st(1,1,:)=c(:,1);
  st(1,2,:)=c(:,2);
  st(2,1,:)=c(:,2);
  st(2,2,:)=c(:,3);
  st(1,3,:)=c(:,4);
  st(3,1,:)=c(:,4);
  st(2,3,:)=c(:,5);
  st(3,2,:)=c(:,5);
  st(3,3,:)=c(:,6);
  
  d = eig3(st)';
  d=real(d);
  % eigenvalues need to be sorted...
  e = zeros(prods,3);
  e(:,1)=min(d,[],2);
  e(:,3)=max(d,[],2);
  e(:,2)=sum(d,2)-e(:,1)-e(:,3);
  d=e;
  
  % New eigenvalues, using J. Weickert's formulas and variants.
  if strcmp(choice,'EED')
      e(:,2:3) = 1.-(1.-alpha)*exp(-(lambda./(e(:,2:3)-e(:,[1,1]))).^m);
      e(:,1)=1;
  elseif strcmp(choice,'cEED')
      e = 1.-(1.-alpha)*exp(-(lambda./e).^m);
  elseif strcmp(choice,'CED')
      e(:,1:2)=alpha+(1.-alpha)*exp(-(lambda./(e(:,[3,3])-e(:,1:2))).^m);
      e(:,3)=alpha;
  elseif  strcmp(choice,'cCED')
      e(:,1:2)=alpha+(1.-alpha)*exp(-((lambda+e(:,1:2))./(e(:,[3,3])-e(:,1:2))).^m);
      e(:,3)=alpha;
  else
      error(['WeichertTensor error : unrecognized '...
          'options.Weickert_choice ' choice ...
          '(EED,cEED,CED,cCED allowed)']);
  end

  
  %interpolating polynomials from old to new eigenvalues.
  dprod = (d(:,1)-d(:,2)).*(d(:,1)-d(:,3)).*(d(:,2)-d(:,3));
  P2 = d(:,2).* e(:,1) - d(:,3).* e(:,1) - d(:,1).*e(:,2) + d(:,3).*e(:,2) + d(:,1).*e(:,3) - d(:,2).*e(:,3);
  P1 = -(d(:,2).^2).*e(:,1) + (d(:,3).^2).* e(:,1) + (d(:,1).^2).*e(:,2) - (d(:,3).^2).*e(:,2) - (d(:,1).^2).*e(:,3) + (d(:,2).^2).* e(:,3);
  P0 = (d(:,2).^2).* d(:,3).* e(:,1) - d(:,2).*( d(:,3).^2).* e(:,1) - (d(:,1).^2).* d(:,3).* e(:,2) + ...
  d(:,1).* (d(:,3).^2) .* e(:,2) + (d(:,1).^2).* d(:,2).* e(:,3) - d(:,1).* (d(:,2).^2).* e(:,3);
  
%    disp(d); disp(e); disp(P0); disp(P1); disp(P2);
  
  %Compute the square of matrices defined by c.
  c2 = zeros([prods,6]);
  c2(:,1)=c(:,1).^2     +c(:,2).^2     +c(:,4).^2;
  c2(:,2)=c(:,1).*c(:,2)+c(:,2).*c(:,3)+c(:,4).*c(:,5);
  c2(:,3)=c(:,2).^2     +c(:,3).^2     +c(:,5).^2;
  c2(:,4)=c(:,1).*c(:,4)+c(:,2).*c(:,5)+c(:,4).*c(:,6);
  c2(:,5)=c(:,2).*c(:,4)+c(:,3).*c(:,5)+c(:,5).*c(:,6);
  c2(:,6)=c(:,4).^2     +c(:,5).^2     +c(:,6).^2;

  
  wt = zeros([prods,6]);
  for i=1:6
    wt(:,i)=(P0*(i==1||i==3||i==6)+P1.*c(:,i)+P2.*c2(:,i))./dprod;
  end
  wt=reshape(wt,[s,6]);  
end

% The alternative construction WeickertTensor_3D_MatlabEig uses Matlab's eig function with a for loop. 
% It should in principle be more robust, in particular when several structure tensor eigenvalues are equal.
% However it is quite slow. (Matlab R2012b)

function tensor = WeickertTensor_3D_MatlabEig(inputTensor,options)

  tensor = inputTensor;
  if(ndims(tensor)~=4 || size(tensor,4)~=6) error('WeickertTensor_3D: invalid tensor format'); end;
  
  lambda = GetOptions(options,'Weickert_lambda');
  alpha = GetOptions(options,'Weickert_alpha',0.1);
  m = GetOptions(options,'Weickert_m',2);
  choice = GetOptions(options,'Weickert_choice','cEED');
  if ~any(strcmp(choice,{'EED','cEED','CED','cCED'}))
      error(['WeichertTensor error : unrecognized '...
          'options.Weickert_choice ' choice ...
          '(EED,cEED,CED,cCED allowed)']);
  end
  
  %preallocating everyone, even the smallest variables...
  s=size(tensor); s=s(1:3); prods = prod(s);
  tensor = reshape(tensor,[prods,6]);
  c = zeros([6,1]); % xx,xy,yy,xz,yz,zz
  mat = zeros([3,3]);  % [xx,xy,xz,xy,yy,yz,xz,yz,zz]
  v=zeros([3,3]);
  d=zeros([3,3]);
  dd=zeros([3,1]);
  
  for i=1:prods % no way to avoid a loop here it seems, due to non-vectorized eig function
      c = tensor(i,:);
      mat=[c(1),c(2),c(4);c(2),c(3),c(5);c(4),c(5),c(6)];
      [v,d]=eig(mat);
      dd=diag(d);
      
      if strcmp(choice,'EED')
          dd(2:3) = 1.-(1.-alpha)*exp(-(lambda./(dd(2:3)-dd(1))).^m);
          dd(1)=1;
      elseif strcmp(choice,'cEED')
          dd = 1.-(1.-alpha)*exp(-(lambda./dd).^m);
      elseif strcmp(choice,'CED')
          dd(1:2)=alpha+(1.-alpha)*exp(-(lambda./(dd(3)-dd(1:2))).^m);
          dd(3)=alpha;
      else % strcmp(choice,'cCED')
          dd(1:2)=alpha+(1.-alpha)*exp(-((lambda+dd(1:2))./(dd(3)-dd(1:2))).^m);
          dd(3)=alpha;
      end
      
    d=diag(dd);
    mat=v*d*v';
    tensor(i,:) = [mat(1,1),mat(1,2),mat(2,2),mat(1,3),mat(2,3),mat(3,3)];
  end

  tensor = reshape(tensor,[s,6]);  
end







