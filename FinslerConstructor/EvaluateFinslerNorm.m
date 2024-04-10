% Copyright Jean-Marie Mirebeau, 2015. jm(dot)mirebeau(at)gmail(dot)com
function result = EvaluateFinslerNorm(u,norm)
    if ndims(u)==2 && size(u,2)==1; u=reshape(u,[1,1,size(u,1)]); end;
    result = sqrt( norm(1)*u(1,:,:).^2 +2*norm(2)*u(1,:,:).*u(2,:,:) +norm(3)*u(2,:,:).^2) - (norm(4).*u(1,:,:)+norm(5).*u(2,:,:));
    result=reshape(result,size(result,2),size(result,3));
end