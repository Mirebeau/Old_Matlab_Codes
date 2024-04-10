% Copyright Jean-Marie Mirebeau, 2015. jm(dot)mirebeau(at)gmail(dot)com
% Transpose the first two levels of a tensor array
% Coincides with usual transposition in dimension 2
function [tx]=Transpose_3D(x)
    d=ndims(x);
    order = 1:d;
    order(1)=2;
    order(2)=1;
    tx=permute(x,order);
end