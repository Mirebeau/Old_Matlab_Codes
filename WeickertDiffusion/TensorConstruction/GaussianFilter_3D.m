% Copyright Jean-Marie Mirebeau, 2015. jm(dot)mirebeau(at)gmail(dot)com

function data = GaussianFilter_3D(data,filter)
    if(~ismatrix(filter)) error('filter must be one dimensional'); end;
    if(size(filter,2)~=1) error('second filter dimension must be singleton'); end;
    data=imfilter(data,         filter,    'symmetric');
    data=imfilter(data,shiftdim(filter,-1),'symmetric');
    data=imfilter(data,shiftdim(filter,-2),'symmetric');
end

function filter = GaussianFilter_3D_Old(w,sigma)
    [x,y,z]=meshgrid(1:w,1:w,1:w);
    o = (w+1)/2.;
    x=x-o;y=y-o; z=z-o;
    filter = exp(-(x.*x+y.*y+z.*z)/(2.*sigma*sigma));
    filter=filter/sum(filter(:));
end