function data = GridCrop(  dims, data )
% Find the intersection of the segments [point, point+offsets] with the boundaries of a box
%  Copyright Jean-Marie Mirebeau, CNRS, University Paris-Sud, 2016

% Input : In substructure data.geometry
%   - offsets? Format : [Dimension,nOffsets] or [Dimension,nPoints,nOffsets]
%   - uniformOffsets (default true)

% Output : In substructure data.geometry
%   - points        Format : [Dimension,nPoints]
%   - multipliers   Format : [1,nPoints]
%   - neighIndices  Format : [1,nPoints]

Dimension = size(dims,2);
nPoints   = prod(dims);
assert(isfield(data.geometry,'offsets'));
offsets = data.geometry.offsets;

if ~(isfield(data.geometry, 'uniformOffsets') && data.geometry.uniformOffsets==false)
    offsets = reshape(offsets,[size(offsets,1),1,size(offsets,2)]);
    offsets = repmat(offsets,[1,nPoints,1]);
end

nOffsets  = size(offsets,3);

assert(size(dims,2)==Dimension);
assert(prod(dims)==nPoints);

multipliers = ones([1,nPoints,nOffsets]);
subscripts = zeros([Dimension,nPoints]);
cOffsets = zeros([1,nPoints,1]);

for iOffset=1:nOffsets

    % Dimension dependent switch is only due to Matlab's syntax.
    switch Dimension
        case 1
            subscripts(1,:) = ind2sub(dims,1:nPoints);
        case 2
            [subscripts(1,:),subscripts(2,:)] = ind2sub(dims,1:nPoints);
        case 3
            [subscripts(1,:),subscripts(2,:),subscripts(3,:)] = ind2sub(dims,1:nPoints);
        otherwise
            disp('Sorry, Grid indices only supports dimension up to 3'); 
            assert(false);
    end
    
    for k=1:Dimension
        cOffsets = offsets(k,:,iOffset);
        pos = cOffsets>0;
        multipliers(1,pos,iOffset)=min( multipliers(1,pos,iOffset), (dims(k)+1-subscripts(k,pos))./cOffsets(pos));
        pos = cOffsets<0;
        multipliers(1,pos,iOffset)=min( multipliers(1,pos,iOffset), -subscripts(k,pos)./cOffsets(pos));
    end
   
end

points = zeros(Dimension,nPoints);
for k=1:Dimension
    mat = 1:dims(k);
    rep = dims; rep(k)=1;
    points(k,:) = reshape(repmat(mat,rep),[1,nPoints]);
end

data.geometry.points = points;
data.geometry.multipliers = multipliers;
data.geometry.neighIndices = GridIndices(dims,offsets);

end

