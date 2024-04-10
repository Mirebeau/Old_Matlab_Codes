% Copyright Jean-Marie Mirebeau, 2015. jm(dot)mirebeau(at)gmail(dot)com
% This function constructs (for each point in a grid) a finsler metric F
% such that : 
% -Forward speed is v. (F(v)=1)
% -Reverse speed ratio is r. (F(-v)=1/r)
% -Side speed ratio is s. (F(R v)=1 where R is the rotation by pi/4)


function finsler = FinslerConstructor(v,r,s)
    squaredNorm = ScalarProduct(v,v);
    V = ScalarVectorDiv(squaredNorm,v);
    a = (1./r + 1.)/2.; a2 = a.^2;
    c = (1./r - 1.)/2.;
    b = sqrt(2)./s + c; b=b.^2 - a2;
    
    m = ScalarVectorProduct(a2,SelfOuterProduct(V)) + ScalarVectorProduct(b, SelfOuterProduct(Perp(V)));
    w = ScalarVectorProduct(c,V);
    
    finsler = cat(1,m,w);
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