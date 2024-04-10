% Copyright Jean-Marie Mirebeau, 2015. jm(dot)mirebeau(at)gmail(dot)com

% Applies to 2-dimensional images of 2-FinslerNorms/3x3-SDP matrices of unit determinant.
% Converts one into the other.

% Format for 2-FinslerNorms : xx,xy,yy,ox,oy, represents pair of 2x2-SDP
% M=[xx,xy;xy,yy] and vector omega = [ox,oy], yielding asymmetric norm
% sqrt{uMu}-omega.u

% Format for 3x3-SDP : xx,xy,yy,xz,yz,zz 

function output = FinslerFromToSDP(input)
    assert(ndims(input)==2);
    assert(size(input,1)==5 || size(input,1)==6);
    s = size(input,2);
    
    if(size(input,1)==5) %FinslerNorms to be converted into 3x3 SDP
        finsler=input;
        sdp = cat(1,finsler,ones([1,s])); %construct representatives
        lambdas = SymmetricMatrix3Determinant(sdp);
%        [minimum,pos]=min(lambdas)
%        sdp(:,pos)
        assert(all(all(lambdas>=0)));
        sdp(6,:)=lambdas; %normalize for unit determinant
        sdp(1,:) = sdp(1,:)./lambdas;
        sdp(2,:) = sdp(2,:)./lambdas;
        sdp(3,:) = sdp(3,:)./lambdas;
        output=sdp; %export
    else
        sdp=input;
        lambdas = sdp(6,:);
        finsler = cat(1,sdp(1,:).*lambdas,sdp(2,:).*lambdas,sdp(3,:).*lambdas,sdp(4,:),sdp(5,:));
        output=finsler;
    end
end