% Returns all pairs of points in p within a given distance d.

% Principle : if |p-q|<=d, then these points lie in a cube
%of side 2r with a center coordinates that are multiples of d.

% Copyright ...

function pairs = PointPairsWithinDistance(p,d)

Dimension = size(p,1);
nPoints = size(p,2);
assert(all(size(p)==[Dimension,nPoints]));

% Make coordinates non-negative
for i=1:Dimension
    p(i,:) = p(i,:)-min(p(i,:));
end

% Round coordinates as multiples of r
q=floor(p/d);

for iOffset = 1:(2^Dimension)
    offset = double(dec2bin(iOffset-1,Dimension)) - double('0');
    for i=1:Dimension
        r(i,:)=floor((q(i,:)+offset(i))/2);
    end
    
    
    % Find points which have identical values of r.
    %(Using sparse matrices ? Sorting ?)
    
    % Find points which are closed than d among these.
end


end