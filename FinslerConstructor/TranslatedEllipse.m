% Copyright Jean-Marie Mirebeau, 2015. jm(dot)mirebeau(at)gmail(dot)com
% This function constructs (for each point in a grid) a finsler metric F
% with prescribed unit ball { x; F(x)<=1 } = { x; (x-W).M.(x-W) <=1 }

function finsler = TranslatedEllipse(M,W)
w=SymmetricMatrixVectorProduct(M,W);
delta = 1.-ScalarProduct(w,W);
m = ScalarVectorDiv(delta.^2, SelfOuterProduct(w)+ScalarVectorProduct(delta,M) );
finsler = cat(1,m, w./delta);
end

% ---------- Basic linear algebra at depth 2 ----------

function w = ScalarVectorProduct(a,v)
    w=v;
    for i=1:size(v,1)
        w(i,:,:) = a.*w(i,:,:);
    end
end

function w = ScalarVectorDiv(a,v)
    w=v;
    for i=1:size(v,1)
        w(i,:,:) = w(i,:,:)./a;
    end
end

function delta = ScalarProduct(v,w)
    delta=sum(v.*w,1);
end

% ----------- Dimension 2 only ------------

function m = SelfOuterProduct(w)
    m = cat(1, w(1,:,:).*w(1,:,:), w(1,:,:).*w(2,:,:), w(2,:,:).*w(2,:,:) );
end

function w = Perp(v)
    w = cat(1, -v(2,:,:), v(1,:,:));
end

function w=SymmetricMatrixVectorProduct(m,v)
    w = cat(1,m(1,:,:).*v(1,:,:)+m(2,:,:).*v(2,:,:), m(2,:,:).*v(1,:,:)+m(3,:,:).*v(2,:,:) );
end