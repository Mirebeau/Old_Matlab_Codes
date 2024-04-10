% Copyright Jean-Marie Mirebeau, 2015. jm(dot)mirebeau(at)gmail(dot)com

function wt = WeickertTensor_2D(tensor,options)
  % Compute the small and large eigenvalues (lambda1,lambda2) of structure tensors
  htr = (tensor(:,:,1)+tensor(:,:,3))/2; %half trace
  det = tensor(:,:,1).*tensor(:,:,3)-tensor(:,:,2).^2;
  sdelta = sqrt(htr.^2 - det);
  lambda1 = htr - sdelta;
  lambda2 = htr + sdelta;
  
%  lambda1
%  lambda2
  
  % Compute the corresponding diffusion tensor eigenvalues (mu1,mu2)
  lambda = GetOptions(options,'Weickert_lambda');
  alpha = GetOptions(options,'Weickert_alpha',0.1);
  m = GetOptions(options,'Weickert_m',2);
  choice = GetOptions(options,'Weickert_choice','cEED');
  s=size(lambda1);
  if strcmp(choice,'EED')
      mu1=ones(s);
      mu2=1.-(1.-alpha)*exp(-(lambda./(lambda2-lambda1)).^m);
  elseif strcmp(choice,'cEED')
      mu1=1.-(1.-alpha)*exp(-(lambda./lambda1).^m);
      mu2=1.-(1.-alpha)*exp(-(lambda./lambda2).^m);
  elseif strcmp(choice,'CED')
      mu1=alpha+(1.-alpha)*exp(-(lambda./(lambda2-lambda1)).^m);
      mu2=alpha;
  elseif strcmp(choice,'cCED')
      mu1=alpha+(1.-alpha)*exp(-((lambda+lambda1)./(lambda2-lambda1)).^m);
      mu2=alpha;
  else error(['WeichertTensor error : unrecognized '...
          'options.Weickert_choice ' choice ...
          '(EED,cEED,CED,cCED allowed)']);
  end
    
%  mu1
%  mu2
  
  % Create the diffusion tensors, with same eigenvectors, new eigenvalues
  wt= ones(size(tensor));
  id=(lambda1==lambda2);
  for i=1:3
      diag= 1-(i==2);
            wt(:,:,i) = (tensor(:,:,i).*(mu1-mu2)-diag*(lambda2.*mu1-lambda1.*mu2))./(lambda1-lambda2+id) ...
          + diag*id.*mu1;
  end
  
%  printf("Weickert ev, tensor. lambda1 : %f, lambda2 : %f, wt %f,%f,%f", lambda1(1,1), lambda2(1,1), wt(1,1,1), wt(1,1,2), wt(1,1,3));
end